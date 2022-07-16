import numpy as np
import h5py
import pprint
import networkx as nx
import matplotlib.pyplot as plt
# from stl import mesh
import pickle
import gzip

import trimesh
import pymeshfix
import glob
import h5py
import pyvista as pv
import itertools
import copy


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
			spine = {'node_dend': id_node_dendrite, 'nodes_terminal': []}
			sub_gs = [g for g in separated_graphs if g.has_node(id_node_dendrite)]
			if len(sub_gs) != 1:
				raise ValueError("multiple/none spine graph(s) was assigned.")
			else:
				sub_g = sub_gs[0]
				ids_node_num_branches = dict( sub_g.degree() )
				spine['nodes_terminal'] = [id for id, num in ids_node_num_branches.items() if num == 1 and id != id_node_dendrite ]
			spines.append(spine)

		return spines


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

		return properties

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
				head.append( {'id': id_head,\
							'node_terminal': node_terminal,\
							'path_from_terminal_to_dend': obtain_path_edges(self.graph, [node_terminal, node_dend]) } )
		return head


	def _obtain_head_properties( self, closed, path ):
		properties = {}
		subgraph   = self.graph.edge_subgraph( path )
		properties['lengths'] = closed.obtain_length_inside( subgraph )
		properties['volumes'] = closed.volume
		properties['areas']   = closed.unclosed_area
		return properties

	def _obtain_neck_properties_straight( self, closed_neck, path ):
		properties = {}
		sub_g = self.graph.edge_subgraph(path)
		properties['lengths'] = closed_neck.obtain_length_inside( sub_g )
		properties['volumes'] = closed_neck.volume
		properties['areas']   = closed_neck.unclosed_area
		return properties

	def _obtain_neck_properties_branched( self, closed_neck, paths ):
		properties = {}
		sub_gs = [ self.graph.edge_subgraph(p) for p in paths ]
		neck_lengths = [closed_neck.obtain_length_inside( sub_g ) for sub_g in sub_gs ]
		ratios       = [l/sum(neck_lengths) for l in neck_lengths]
		properties['lengths'] = neck_lengths
		properties['volumes'] = [closed_neck.volume*ratio for ratio in ratios]
		properties['areas']   = [closed_neck.unclosed_area*ratio for ratio in ratios]
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


	def create( self ):
		# Find the spine heads associated with terminal nodes, then it is associated with a spine-neck.
		for spine in self.spines:
			nodes_terminal = spine['nodes_terminal']
			node_dend      = spine['node_dend']
			spine['head']  = self._associate_ids_head_and_nodes_terminal(node_dend, nodes_terminal)
			spine['neck']  = [ {'id': self.adjacents.obtain_a_neck_connecting_to_head( head['id'] ) } for head in spine['head'] ]

		# Obtain the compartment of each spine head.
		for spine in self.spines:
			for head in spine['head']:
				closed_head = ClosedMesh( self.vertices, self.head_fs[head['id']] )
				head.update( self._obtain_head_properties( closed_head,  head['path_from_terminal_to_dend']) )

		# Aggregate the neck that has multiple heads
		for spine in self.spines:
			ids_neck  = [n['id'] for n in spine['neck']]
			if len(ids_neck) != len(set(ids_neck)):
				unique_necks, heads_for_unique_necks = self._aggregate_neck_with_mutiple_heads( ids_neck, spine['head'] )
				spine['neck'] = [{'id': id } for id in unique_necks ]
				spine['head'] = heads_for_unique_necks

		# Obtain the compartments of each spine neck.
		for spine in self.spines:
			for neck, heads in zip( spine['neck'], spine['head'] ):
				closed_neck = ClosedMesh( self.vertices, self.neck_fs[ neck['id'] ] )
				nodes       = closed_neck.obtain_nodes_inside( self.graph )
				if isinstance(heads, dict): # Single head
					neck_properties = self._obtain_neck_properties_straight( closed_neck, heads['path_from_terminal_to_dend'] ) 
					neck.update( neck_properties )
				elif len(nodes) == 0:
					ids_head = [ h['id'] for h in heads ]
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
					ids_head = [ h['id'] for h in heads ]
					print('heads {} are the member of neck_id {} that has the nodes {}'.format(ids_head, neck['id'], nodes))
					print('Neck shape is complicated beyond the current implementation (len(nodes) > 1). Ignored.')





#
class AdjacentFaces():
	def __init__(self, vertices, head_fs, neck_fs):
		self.heads_boundary_vertices = [obtain_boundary_vertices(vertices, f) for f in head_fs]
		self.necks_boundary_vertices = [obtain_boundary_vertices(vertices, f) for f in neck_fs]

	def obtain_a_neck_connecting_to_head(self, id_head):
		id_neck = self.obtain_necks_connecting_to_head( id_head )
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

	def obtain_necks_connecting_to_head(self, id_head):
		head_vertices = self.heads_boundary_vertices[id_head]
		adjacents = []
		for i, neck_vertices in enumerate(self.necks_boundary_vertices):
			if neck_vertices & head_vertices != set():
				adjacents.append(i)
		
		return adjacents

def obtain_boundary_vertices(vertices, faces):

	if faces == []:
		return set()

	mesh = trimesh.Trimesh(vertices, np.array(faces), process=False)
	#mesh.remove_degenerate_faces()
	mesh.remove_duplicate_faces()
	unique_edges = mesh.edges[trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)]
	boundary_vertices = set(unique_edges.flatten())
	return boundary_vertices


def load_paint_subtraction(filenames_for_paint_subtraction, faces):

	# Subtract reference meshes
	subtracted_face_Flag = [True]*faces.shape[0]
	for filename in filenames_for_paint_subtraction:
		with open(filename, 'rb') as file:
			data = pickle.load(file)
		data = data['painted']
		unzipped_tri = gzip.decompress(data)
		for i in range( faces.shape[0] ) :
			if (unzipped_tri[i*3:i*3+3] == b'\x01\x01\x01') :
				subtracted_face_Flag[i] = False

	# Get faces (Ids of three vertices)
	subtracted_face = []
	for i in range( faces.shape[0] ) :
		if subtracted_face_Flag[i] == True:
			subtracted_face.append(faces[i,:])
	subtracted_face = np.array(subtracted_face)

	# Obtain separated faces
	subtracted_faces = _separate_unconnected_faces(subtracted_face, num_faces_to_remove_as_fragumented_face = 100)
	return subtracted_faces


##
class ClosedMesh():
	def __init__(self, vertices, faces):

		part_mesh = pymeshfix.MeshFix(vertices, faces)
		part_mesh.repair()
		self.mesh = trimesh.Trimesh(part_mesh.v, part_mesh.f)
		self.mesh.merge_vertices()
		self.mesh.remove_degenerate_faces()
		self.mesh.remove_duplicate_faces()
		self.volume = self.mesh.volume
		
		mesh0 = trimesh.Trimesh(vertices, faces)
		self.unclosed_area = mesh0.area

	def obtain_crossing_edges(self, graph):
		crossing_edges = []
		for id_edge in graph.edges.keys():
			loc0 = graph.nodes[id_edge[0]]['loc'].reshape([1,-1])
			loc1 = graph.nodes[id_edge[1]]['loc'].reshape([1,-1])
			in0  = self.mesh.contains(loc0)
			in1  = self.mesh.contains(loc1)
			if in0 != in1:
				crossing_edges.append(id_edge)
		return crossing_edges

	def obtain_length_inside(self, graph):
		return sum( [self.obtain_inside_length_edge( edge ) for edge in graph.edges.values()] )

	def obtain_nodes_inside(self, graph):
		nodes_inside = []
		for id_node in graph.nodes.keys():
			loc0 = graph.nodes[id_node]['loc'].reshape([1,-1])
			if self.mesh.contains(loc0):
				nodes_inside.append(id_node)
		return nodes_inside

	def obtain_inside_length_edge(self, edge):
		path    = edge['path']
		# points ((n, 3) float)
		# return (n, ) bool
		enclosed = self.mesh.contains(path)
		sub_path = path[enclosed == True,:]
		
		diff_sub_path   = np.diff(sub_path, axis=0)
		length_sub_path = np.sum(np.linalg.norm(diff_sub_path, axis=1))
		#print('length_sub_path ', length_sub_path)
		return length_sub_path


'''
def obtain_necks_connecting_to_each_heads(vertices, head_fs, neck_fs):

		heads_boundary_vertices = [obtain_boundary_vertices(vertices, f) for f in head_fs]
		necks_boundary_vertices = [obtain_boundary_vertices(vertices, f) for f in neck_fs]
		neck_connecting_to_neck = []
		for i, head_vertices in enumerate(heads_boundary_vertices):
		for i, neck_vertices in enumerate(necks_boundary_vertices):
			adjacents = []
			for j, head_vertices in enumerate(heads_boundary_vertices):
				#print('len(ref_boundary_vertices) ', len(ref_boundary_vertices))
				if neck_vertices & head_vertices != set():
					adjacents.append(j)
			heads_bind_to_neck.append(adjacents)
		print('heads_bind_to_neck ', heads_bind_to_neck )
'''


def obtain_mesh_barycenter(vertices, faces):
	mesh       = trimesh.Trimesh( vertices, faces )
	vertices   = mesh.vertices
	barycenter = np.mean( vertices, axis=0 )
	#print('barycenter.shape ', barycenter.shape )
	return barycenter


def load_stl(filename_mesh):
	with open(filename_mesh, 'rb') as f:
		mesh = trimesh.exchange.stl.load_stl(f) # trimesh.exchange.stl.load_stl_binary
		vertices = mesh['vertices']
		faces    = mesh['faces']

	mesh     = trimesh.Trimesh(vertices, faces)
	vertices = mesh.vertices
	faces    = mesh.faces
	return vertices, faces


'''
def load_stl_old(filename_mesh):
	m = mesh.Mesh.from_file(filename_mesh)
	shape = m.points.shape
	vertices = m.points.reshape(-1, 3)
	faces    = np.arange(vertices.shape[0]).reshape(-1, 3)
	return vertices, faces
'''


def _separate_unconnected_faces(faces, num_faces_to_remove_as_fragumented_face = 20 ):

	# Obtain separated faces
	if type(faces).__module__ == "numpy":
		faces = faces.tolist()
	if faces == []:
		return []
	# print('faces ', faces)
	adjacency = trimesh.graph.face_adjacency(faces=faces, mesh=None, return_edges=False)
	graph = nx.Graph()
	graph.add_edges_from(adjacency)
	ids_face = list( nx.connected_components(graph) )

	# Separate unconnected faces
	separated_faces    = []
	boundary_vertices = []
	faces = np.array( faces )
	for id in ids_face:
		#print('In load_paint_subtraction, len(face_id): ', len(id))
		if len(id) < num_faces_to_remove_as_fragumented_face:
			continue
		# print("list(id) ", list(id))
		# print("id ", id)
		separated_faces.append(faces[list(id),:] )
	return separated_faces


def _obtain_paint_from_file(filename_paint, faces):
	with open(filename_paint, 'rb') as file:
		data = pickle.load(file)
	data = data['painted']

	unzipped_tri = gzip.decompress(data)
	sub_face = []
	for i in range( faces.shape[0] ) :
		if (unzipped_tri[i*3:i*3+3] == b'\x01\x01\x01') :
			sub_face.append(faces[i,:])
	return np.array(sub_face)


def load_paints(filenames_paint, faces):
	multiple_faces = []
	for filename_paint in filenames_paint:
		sub_face = _obtain_paint_from_file(filename_paint, faces)
		multiple_faces.extend( _separate_unconnected_faces(sub_face) )
	return multiple_faces


def load_paint(filename_paint, vertices, faces):
	sub_face = _obtain_paint_from_file(filename_paint, faces)
	if sub_face == []:
		vclean, fclean, fcenter, farea = None, None, None, None
		return vclean, fclean, fcenter, farea

	mesh = trimesh.Trimesh(vertices, np.array(sub_face))
	mesh.merge_vertices()
	mesh.remove_degenerate_faces()
	mesh.remove_duplicate_faces()
	vclean  = np.array(mesh.vertices)
	fclean  = np.array(mesh.faces)
	fcenter = np.array(mesh.triangles_center)
	farea   = np.array(mesh.area_faces)
	return vclean, fclean, fcenter, farea


def fill_hole(vertices, faces):
	part_mesh = pymeshfix.MeshFix(vertices, faces)
	part_mesh.repair()

	closed_v = part_mesh.v # numpy np.float array
	closed_f = part_mesh.f # numpy np.int32 array
	
	mesh = trimesh.Trimesh(closed_v, closed_f)
	mesh.merge_vertices()
	mesh.remove_degenerate_faces()
	mesh.remove_duplicate_faces()
	vclean  = np.array(mesh.vertices)
	fclean  = np.array(mesh.faces)
	fcenter = np.array(mesh.triangles_center)
	farea   = np.array(mesh.area_faces)
	volume  = mesh.volume
	return vclean, fclean, fcenter, farea, volume


class CreateGraph():
	"""Create a graph from hdf5.

	Args:
		filename_hdf5 (str): Target Hdf5 filename

	Returns:
		(pyLD.CreateGraph): CreateGraph object that has the follwing instances:
		- graph (float): Graph in a networkx object.
		- node_src_id_dst_id (list): List contains two integers that denotes the ids of start-and-end nodes.
	"""
	def __init__(self, filename_hdf5):
		self._load_skeleton_file(filename_hdf5)
		self.ids_edge = np.unique( self.org_ids_edge )
		edges = self.obtain_edge()
		nodes, self.ids_node = self.obtain_node(edges)
		edges = self.associate_edge_with_node(nodes, edges)
		self.graph = self.make_graph(nodes, edges)
		self.node_src_id_dst_id = self.find_dendrite()

	def obtain_edge(self):
		edges = {}
		for id in self.ids_edge:
			# Obtain line
			org_targ_edges = self.org_edges[self.org_ids_edge == id, :]
			pos = self.org_vertices[org_targ_edges[:,0], :]
			pos = np.vstack( [pos, self.org_vertices[org_targ_edges[-1,1], :]])
			
			tangents = self.org_tangents[org_targ_edges[:,0], :]
			tangents = np.vstack( [tangents, self.org_tangents[org_targ_edges[-1,1], :]])
			
			radiuses = self.org_radiuses[org_targ_edges[:,0]]
			radiuses = np.hstack( [radiuses, self.org_radiuses[org_targ_edges[-1,1]]])
			
			# Obtain length
			pos_diff   = np.diff(pos, axis=0)
			pos_length = np.sum(np.linalg.norm(pos_diff, axis=1))
			# Save properties
			edges[id] = {'path': pos,\
				'tangents': tangents,\
				'radiuses': radiuses,\
				'start': pos[0,:],\
				'end': pos[-1,:],\
				'length': pos_length,\
				'nodes': []}
		return edges


	def _obtain_edge_points(self, edges):
		#
		edge_points = np.empty([0,5])
		for id in self.ids_edge:
			start_p = np.hstack( [edges[id]['start'], 0, id] )
			end_p   = np.hstack( [edges[id]['end']  , -1, id ] )
			points = np.vstack( [ start_p, end_p ] )
			edge_points = np.vstack([edge_points, points])
		return edge_points
		# edge_points = 
		#     np.array([[x_start, y_start, z_start, 0, id1],
		#               [x_end, y_end, z_end, -1, id1],
		#               [,,,id2],...])


	def obtain_node(self, edges):
		# Obtain edge points
		edge_points = self._obtain_edge_points(edges)
		# Find shared start-end points.
		connections = []
		for id in range(edge_points.shape[0]):
			current_edge_points = edge_points[id,:]
			diff     = edge_points[:,:3] - current_edge_points[:3]
			diff_len = np.linalg.norm(diff, axis=1)
			tmp_locs = np.hstack( current_edge_points )
			tmp_locs = tmp_locs[np.newaxis,:]
			for id_len in np.argsort(diff_len)[1:]: # Without its own point
				if (diff_len[id_len] < 0.07):       # if they are suffficiently close
					if(int( edge_points[id_len,4] ) not in tmp_locs[:,4].tolist() ): # No double edges
						tmp_locs = np.vstack( [ tmp_locs, edge_points[id_len,:] ]  )
				else:
					connection = tmp_locs[np.argsort(tmp_locs[:, 4])]
					connections.append( connection )
					break
		
		# Remove multiplicated connections and obain nodes
		# print('connections ', connections)
		id_node = 0
		nodes    = {}
		nodes[id_node] = {'raw': connections[0], \
			'edges': connections[0][:,4].astype('int').tolist(), \
			'edges_start_0_end_m1': connections[0][:,3].astype('int').tolist(), \
			'loc': np.mean(connections[0][:,:3], axis=0)}
		for con in connections:
			flag = 0
			for unique_con in nodes.values():
				if (con.shape[0] == unique_con['raw'].shape[0]) and (np.all(con - unique_con['raw'] < 1.0e-1) ):
					flag = 1
					break
			if (flag == 0):
				nodes[id_node] = {'raw': con, \
					'edges': con[:,4].astype('int').tolist(), \
					'edges_start_0_end_m1': con[:,3].astype('int').tolist(), \
					'loc': np.mean(con[:,:3], axis=0)}
				id_node += 1
		ids_node = list(range(id_node))
		print('len(connections) ', len(connections))
		print('len(nodes)       ', len(nodes))
		return nodes, ids_node


	def associate_edge_with_node(self, nodes, edges):
		# Associate edges with the nodes
		for id_node in self.ids_node:
			for id_edge in nodes[id_node]['edges']:
				edges[id_edge]['nodes'].append(id_node)
		return edges


	def make_graph(self, nodes, edges):
		graph = nx.Graph()
		for id_node, node in nodes.items():
			graph.add_node(id_node, \
				raw                  = node['raw'], \
				edges                = node['edges'], \
				edges_start_0_end_m1 = {i:v for i,v in zip(node['edges'], node['edges_start_0_end_m1'])}, \
				loc                  = node['loc']
				)
		
		for id_edge, edge in edges.items():
			if ( len(edge['nodes']) == 2 ):
				graph.add_edge(edge['nodes'][0], edge['nodes'][1], \
					id       = id_edge, \
					len      = edge['length'], \
					path     = edge['path'], \
					tangents = edge['tangents'], \
					radiuses = edge['radiuses']
					)
		return graph
		# http://47.112.232.56/a/stackoverflow/en/629b5e10652087181359e77a.html


	def find_dendrite(self):
		id_src_id_dst_dist = []
		for id_source, ddistances in nx.all_pairs_dijkstra_path_length(self.graph,  weight='len'):
			max_target_id_dist   = max( ddistances.items(), key=(lambda x: x[1]))
			id_src_id_dst_dist.append( [id_source, max_target_id_dist[0], max_target_id_dist[1]] )
		max_id_src_id_dst_dist = max( id_src_id_dst_dist, key=(lambda x: x[2]) )
		#print('src_id_dst_idt: ', max_id_src_id_dst_dist[:2], ' , distance: ', max_id_src_id_dst_dist[2])
		return max_id_src_id_dst_dist[:2]


	def _load_skeleton_file(self, filename_hdf5):
		with h5py.File( filename_hdf5 ,'r') as f:
			self.org_ids_edge = f['ids_edge'][()]
			self.org_edges    = f['edges'][()]
			self.org_radiuses = f['radiuses'][()]
			#self.org_lengths  = f['lengths'][()]
			self.org_vertices = f['vertices'][()]
			self.org_tangents = f['tangents'][()]


def obtain_nodes_location(graph, ids_node):
	locations = []
	for id_node in ids_node:
		loc = graph.nodes[id_node]['loc']
		locations.append( loc )
	locations = np.array(locations)
	return locations


#def obtain_ids_graph_edge_dendrite(graph, src_id_dst_id):
#	return list(nx.all_simple_edge_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]

#def obtain_ids_graph_node_dendrite(graph, src_id_dst_id):
#	return list(nx.all_simple_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]

def obtain_path_nodes(graph, src_id_dst_id):
	return list(nx.all_simple_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]

def obtain_path_edges(graph, src_id_dst_id):
	return list(nx.all_simple_edge_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]


def obtain_distance_between_edges_dendrite(graph, src_id_dst_id):
	# Edge lines
	lengths_edges_dendrite = []
	for id_graph_edge in obtain_path_edges(graph, src_id_dst_id):
		loc0  = graph.nodes[id_graph_edge[0]]['loc']
		loc1  = graph.nodes[id_graph_edge[1]]['loc']
		lengths_edges_dendrite.append( np.linalg.norm(loc0-loc1) )
	return lengths_edges_dendrite


def obtain_ids_edge_dendrite(graph, src_id_dst_id):
	# Edge lines
	ids_edge_dendrite     = []
	for id_graph_edge in obtain_path_edges(graph, src_id_dst_id):
		id    = graph.edges[id_graph_edge]['id']
		ids_edge_dendrite.append( id )
	return ids_edge_dendrite


def obtain_locs_tangents_from_dendrite(graph, src_id_dst_id):
	#
	ids_graph_node_dendrite = obtain_path_nodes(graph, src_id_dst_id)

	# Tangents
	locations = obtain_nodes_location(graph, ids_graph_node_dendrite)
	loc2  = np.delete( locations, [-1, -2], 0)
	loc1  = np.delete( locations, [0, 1], 0)
	tangents  = loc2-loc1
	n_conv = 4
	k = np.ones(n_conv)/n_conv
	for i in range(3):
		tangents[:,i] = np.convolve(tangents[:,i], k, mode='same')
	#print(tangents.shape[0])
	for j in range(tangents.shape[0]):
		tangents[j,:] = tangents[j,:]/np.linalg.norm( tangents[j,:] )
	return locations, tangents

def obtain_ids_of_terminal_nodes_from_dendrite(graph, src_id_dst_id):
	#
	# Renumber terminal nodes along the dendrite
	ids_graph_node_dendrite = obtain_path_nodes(graph, src_id_dst_id)
	ids_node_num_branches = dict( graph.degree() )
	ids_node_terminal     = [i for i, num in ids_node_num_branches.items() if num == 1 and not i in ids_graph_node_dendrite ]
	distances             = dict(nx.all_pairs_shortest_path_length(graph))
	dendritic_node_has    = {i: [] for i in ids_graph_node_dendrite }
	for id_node_terminal in ids_node_terminal:
		dists = {id_node_dend: distances[id_node_terminal][id_node_dend] for id_node_dend in ids_graph_node_dendrite}
		dendritic_node_has[min(dists, key=dists.get)].append(id_node_terminal)
	#print('dendritic_node_has ', dendritic_node_has)
	
	return dendritic_node_has

# '''


# '''


class NodeBinaryTree:
	def __init__(self, data):
		self.left = None
		self.right = None
		self.data = data
	def preorder(node):
		if node:
			print(node.data)
			preorder(node.left)
			preorder(node.right)


	'''
	# Renumber end-nodes along the dendrite
	nodes_num_branches = dict( graph.degree() )
	
		[for id_dend in ids_graph_node_dendrite]
	nodes_end     = [i for i, num in nodes_num_branches.items() if num == 1 and i not in src_id_dst_id] #  and i not in src_id_dst_id
	nodes_end_len = [len(list(nx.all_simple_paths(graph, src_id_dst_id[0], node_end))[0]) for node_end in nodes_end]
	ids_len       = sorted(range(len(nodes_end_len)), key=lambda k: nodes_end_len[k])
	print('ids_len ', ids_len)
	nodes_end     = [nodes_end[i] for i in ids_len]
	'''


if __name__ == '__main__':
	filename_hdf5 = 'data2/0000000001.hdf5'
	a = CreateGraph(filename_hdf5)
	nx.draw(a.graph, with_labels=True)
	plt.show()

