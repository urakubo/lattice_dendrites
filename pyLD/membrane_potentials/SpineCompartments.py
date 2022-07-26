import numpy as np
import networkx as nx
import trimesh
import pymeshfix
import itertools
import copy

from CloseMesh  import CloseMesh
from Utils import obtain_path_nodes, obtain_path_edges


class SpineCompartments():

	def __init__( self, graph, node_src_id_dst_id, vertices, head_fs, neck_fs, barycenters_head ):
		# Preparation
		self.graph        = graph
		self.spines       = self._obtain_putative_spine_terminals_along_dendrite( node_src_id_dst_id )
		self.adjacents    = AdjacentFaces( vertices, head_fs, neck_fs )
		self.vertices     = vertices
		self.head_fs      = head_fs
		self.neck_fs      = neck_fs
		self.barycenters_head = barycenters_head


	def create( self ):
		# Find the spine heads associated with terminal nodes, then each head is associated with a spine-neck.
		# Terminal-based FOR
		for spine in self.spines:
			nodes_terminal = spine['nodes_terminal']
			node_dend      = spine['node_dend']
			spine['head']  = self._associate_ids_head_and_nodes_terminal(node_dend, nodes_terminal)
			spine['neck']  = [ {'id_neck': self.adjacents.obtain_a_neck_connecting_to_head( head['id_head'] ) } for head in spine['head'] ]

		# Obtain the compartment of each spine head.
		# Head-based FOR
		for spine in self.spines:
			for head in spine['head']:
				closed_head = CloseMesh( self.vertices, self.head_fs[head['id_head']] )
				head.update( self._obtain_head_properties( closed_head,  head['path_from_terminal_to_dend']) )

		# Aggregate the neck that has multiple heads
		for spine in self.spines:
			ids_neck  = [n['id_neck'] for n in spine['neck']]
			if len(ids_neck) != len(set(ids_neck)):
				unique_necks, heads_for_unique_necks = self._aggregate_neck_with_mutiple_heads( ids_neck, spine['head'] )
				spine['neck'] = [{'id_neck': id } for id in unique_necks ]
				spine['head'] = heads_for_unique_necks

		# Obtain the compartments of each spine neck.
		# Neck-based FOR
		for spine in self.spines:
			for neck, heads in zip( spine['neck'], spine['head'] ):
				closed_neck = CloseMesh( self.vertices, self.neck_fs[ neck['id_neck'] ] )
				nodes       = closed_neck.obtain_nodes_inside( self.graph )
				if isinstance(heads, dict): # Single head
					neck_properties = self._obtain_neck_properties_straight( closed_neck, heads['path_from_terminal_to_dend'] ) 
					neck.update( neck_properties )
				elif len(nodes) == 0:
					ids_head = [ h['id_head'] for h in heads ]
					paths    = [ h['path_from_terminal_to_dend'] for h in heads ]
					print('Nearly separated branching neck (Neck id: {}, Head ids: {}. Spine neck is divided depending on its length.'.format(neck['id'], ids_head))
					neck.update( self._obtain_neck_properties_branched( closed_neck, paths ) )
				elif len(nodes) == 1:
					print('Branching spines with the node', nodes[0])
					paths   = [ h['path_from_terminal_to_dend'] for h in heads ]
					paths_Y = self._obtain_paths_compartment_from_Y_shaped_paths( paths )
					neck.update( self._obtain_neck_properties_branched( closed_neck, paths_Y ) )
					neck['node'] = nodes[0]
				else:
					ids_head = [ h['id_head'] for h in heads ]
					print('heads {} are the member of neck_id {} that has the nodes {}'.format(ids_head, neck['id_neck'], nodes))
					print('Neck shape is complicated beyond the current implementation (len(nodes) > 1). Ignored.')


	def _obtain_putative_spine_terminals_along_dendrite(self, src_id_dst_id):
		'''Obtain graphs from each node of dendritic shaft.

		Args:
			src_id_dst_id (list[int]): Source and destination ids of nodes of the dendritc shaft
		
		Returns:
			spines (list[dict]): 'node_dend': dendritic node, 'nodes_terminal': terminal notes without the dendritic node
		'''

		ids_graph_edge_dendrite = obtain_path_edges(self.graph, src_id_dst_id)
		graph_separated = copy.deepcopy(self.graph)
		graph_separated.remove_edges_from(ids_graph_edge_dendrite)
		separated_graphs = [graph_separated.subgraph(c) for c in nx.connected_components(graph_separated)]

		ids_graph_node_dendrite = obtain_path_nodes(self.graph, src_id_dst_id)
		spines = []
		for i, id_node_dendrite in enumerate(ids_graph_node_dendrite):
			spine  = {'node_dend': id_node_dendrite, 'nodes_terminal': []}
			sub_gs = [g for g in separated_graphs if g.has_node(id_node_dendrite)]
			if len(sub_gs) != 1:
				raise ValueError("multiple/none spine graph(s) was assigned.")
			else:
				sub_g = sub_gs[0]
				ids_node_num_branches = dict( sub_g.degree() )
				spine['nodes_terminal'] = [id for id, num in ids_node_num_branches.items() if num == 1 and id != id_node_dendrite ]
			spines.append(spine)
		return spines


	def _obtain_neck_properties_branched( self, closed_neck, paths ):
		properties = {}
		sub_gs = [ self.graph.edge_subgraph(p) for p in paths ]
		neck_lengths = closed_neck.obtain_lengths_inside( sub_gs )
		faces_for_branch, fareas_for_branch = closed_neck.assign_mesh_for_each_graph( sub_gs ) ### face assignments
		ratios       = [l/sum(neck_lengths) for l in neck_lengths] ### length ratios
		properties['lengths_for_branch'] = neck_lengths
		properties['volumes_for_branch'] = [closed_neck.volume*ratio for ratio in ratios] ### 円筒だと思うと Volume = S^2 / (4*pi*L)
		properties['faces_for_branch']   = faces_for_branch
		properties['areas_for_branch']   = fareas_for_branch
		return properties


	def _aggregate_neck_with_mutiple_heads( self, ids_neck, heads ):
		ids_unique_necks, inverse, counts = np.unique( ids_neck, return_inverse=True, return_counts=True )
		inverse = inverse.tolist()
		heads_for_unique_necks = [[] for _ in range(ids_unique_necks.shape[0])]
		for id, head in zip(inverse, heads):
			if counts[id] == 1:
				heads_for_unique_necks[id] = head
			else:
				heads_for_unique_necks[id].append(head)
		return ids_unique_necks, heads_for_unique_necks


	def _associate_ids_head_and_nodes_terminal( self, node_dend, nodes_terminal ):
		head     = []
		ids_head = []
		for node_terminal in nodes_terminal:
			loc_terminal = self.graph.nodes[node_terminal]['loc']
			dists  = np.linalg.norm(self.barycenters_head-loc_terminal, axis=1)
			id_head = np.argmin(dists)
			if id_head in ids_head:
				print('Warning: the selected spine head is again selected. The latter is ignored.')
			else:
				ids_head.append(id_head)
				head.append( {'id_head': id_head,\
							'node_terminal': node_terminal,\
							'path_from_terminal_to_dend': obtain_path_edges(self.graph, [node_terminal, node_dend]) } )
		return head


	def _obtain_head_properties( self, closed, path ):
		properties = {}
		subgraph   = self.graph.edge_subgraph( path )
		properties['length'] = closed.obtain_length_inside( subgraph )
		properties['volume'] = closed.volume
		properties['area']   = closed.unclosed_area
		return properties


	def _obtain_neck_properties_straight( self, closed_neck, path ):
		properties = {}
		sub_g = self.graph.edge_subgraph(path)
		properties['length'] = closed_neck.obtain_length_inside( sub_g )
		properties['volume'] = closed_neck.volume
		properties['area']   = closed_neck.unclosed_area
		return properties


	def _obtain_paths_compartment_from_Y_shaped_paths( self, paths ):
		"""Obtain paths for each head at first, then obtain a shared path.

		Args:
			paths (list[path1, path2,..., path(n)]): List of edged-based paths from_terminals_to_a_dend

		Returns:
			(list[path_head1, path_head2, ..., path_head(n), path_shared]): Separated paths.
		"""
		connected_path         = sum(paths,[])
		total_path             = set(connected_path)
		one_time_used_paths    = set([x for x in total_path if connected_path.count(x) == 1])
		path_specific_paths    = [list( set(p) & one_time_used_paths ) for p in paths]
		any_time_used_paths    = [x for x in total_path if connected_path.count(x) == len(paths)]
		arranged_paths = path_specific_paths + [any_time_used_paths]

		return arranged_paths


#
class AdjacentFaces():
	def __init__(self, vertices, head_fs, neck_fs):
		self.heads_boundary_vertices = [self._obtain_boundary_vertices(vertices, f) for f in head_fs]
		self.necks_boundary_vertices = [self._obtain_boundary_vertices(vertices, f) for f in neck_fs]

	def obtain_a_neck_connecting_to_head(self, id_head):
		id_neck = self._obtain_necks_connecting_to_head( id_head )
		if len(id_neck) >= 2:
			print('Warning: more than two necks {} are connected to a head {}. The first one is selected for connection.'.format(id_neck, id_head))
			id_neck = id_neck[0]
		elif len(id_neck) == 0:
			print('id_head ', id_head)
			print('Warning: no spine neck is connected to a head.')
			id_neck = None
		else :
			id_neck = id_neck[0]

		return id_neck

	def _obtain_necks_connecting_to_head(self, id_head):
		head_vertices = self.heads_boundary_vertices[id_head]
		adjacents = []
		for i, neck_vertices in enumerate(self.necks_boundary_vertices):
			if neck_vertices & head_vertices != set():
				adjacents.append(i)
		
		return adjacents

	def _obtain_boundary_vertices(self, vertices, faces):

		if faces == []:
			return set()
		mesh = trimesh.Trimesh(vertices, np.array(faces), process=False)
		#mesh.remove_degenerate_faces()
		mesh.remove_duplicate_faces()
		unique_edges = mesh.edges[trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)]
		boundary_vertices = set(unique_edges.flatten())
		return boundary_vertices

