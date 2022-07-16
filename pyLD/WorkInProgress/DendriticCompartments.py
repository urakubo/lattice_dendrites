import numpy as np
import networkx as nx
import trimesh
import itertools

from CloseMesh  import CloseMesh

class DendriticCompartments():
	def __init__( self, graph, vertices_org, faces_org ):
		# Preparation
		self.graph  = graph
		self.v_org  = vertices_org
		self.f_org  = faces_org
		m           = trimesh.Trimesh(self.v_org, self.f_org)
		self.fcenters_org = np.array(m.triangles_center)
		self.fareas_org   = np.array(m.area_faces)
		
		closed_d = CloseMesh(self.v_org, self.f_org)
		self.v_fill        = closed_d.obtain_vertices()
		self.f_fill        = closed_d.obtain_faces()
		self.fcenters_fill = closed_d.obtain_fcenters()

	def update_nodes( self, src_id_dst_id ):
		self.nodes     = list(nx.all_simple_paths(self.graph, src_id_dst_id[0], src_id_dst_id[1]))[0]
		self.locations = self._obtain_nodes_location(self.nodes)
		self.tangents  = self._obtain_nodes_tangents(self.nodes)
		self.distances = self._obtain_distance_between_edges_dendrite(self.nodes)

	def _obtain_distance_between_edges_dendrite(self, nodes):
		# location
		locs = []
		for node in nodes:
			locs.append(self.graph.nodes[node]['loc'])
		locs  = np.array(locs)
		diffs = np.diff(locs, axis=0)
		dists = np.linalg.norm(locs, axis=1)
		return dists

	def _obtain_nodes_location(self, nodes):
		locations = [self.graph.nodes[node]['loc'] for node in nodes]
		locations = np.array(locations)
		return locations

	def _obtain_nodes_tangents(self, nodes):
		# Tangents
		locations = self._obtain_nodes_location(nodes)
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
		return tangents

	def _split_faces(self, fcenters):
		locations_ = np.delete( self.locations, [0, -1], 0)
		ids_f_part = [ np.dot(fcenters-loc, tang) > 0 for tang, loc in zip(self.tangents, locations_)]

		i_part_f = []
		for i in range(len(self.tangents)+1):
			if i == 0:
				i_part_f.append( ids_f_part[i] )
			elif i < len(self.tangents):
				i_part_f.append( ids_f_part[i] * np.logical_not(ids_f_part[i-1]) )
			elif i == len(self.tangents):
				i_part_f.append( np.logical_not(ids_f_part[i-1]) )
			else:
				print('Error.')
				break
		return i_part_f

	def create(self):
		# Obtain faces
		face_ids = self._split_faces( self.fcenters_org  )
		areas    = [ np.sum( self.fareas_org[ids] ) for ids in face_ids ]

		face_ids = self._split_faces( self.fcenters_fill )
		closed   = [CloseMesh(self.v_fill, self.f_fill[ids]) for ids in face_ids ]
		vertices = [c.obtain_vertices() for c in closed ]
		faces    = [c.obtain_faces()    for c in closed ]
		volumes  = [c.volume            for c in closed ]

		return vertices, faces, volumes, areas


