import numpy as np
import h5py
import pprint
import networkx as nx
import matplotlib.pyplot as plt
from stl import mesh
import pickle
import gzip

import trimesh
import pymeshfix
import glob
import h5py
import pyvista as pv
import pymeshfix


def load_paint_subtraction(filenames_paint_subtraction, vertices, faces):

	sub_face_Flag = [True]*faces.shape[0]

	for filename in filenames_paint_subtraction:
		with open(filename, 'rb') as file:
			data = pickle.load(file)
		data = data['painted']
		unzipped_tri = gzip.decompress(data)
		for i in range( faces.shape[0] ) :
			if (unzipped_tri[i*3:i*3+3] == b'\x01\x01\x01') :
				sub_face_Flag[i] = False

	sub_face = []
	for i in range( faces.shape[0] ) :
		if sub_face_Flag[i] == True:
			sub_face.append(faces[i,:])
	
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
		


def load_stl(filename_mesh):
	m = mesh.Mesh.from_file(filename_mesh)
	shape = m.points.shape
	vertices = m.points.reshape(-1, 3)
	faces    = np.arange(vertices.shape[0]).reshape(-1, 3)
	return vertices, faces

def load_paint(filename_paint, vertices, faces):
	with open(filename_paint, 'rb') as file:
		data = pickle.load(file)
	data = data['painted']

	unzipped_tri = gzip.decompress(data)
	sub_face = []
	for i in range( faces.shape[0] ) :
		if (unzipped_tri[i*3:i*3+3] == b'\x01\x01\x01') :
			sub_face.append(faces[i,:])

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
	"""
	Create a graph from hdf5.

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


def obtain_edges_dendrite(graph, src_id_dst_id):
	#
	ids_graph_edge_dendrite = \
		list(nx.all_simple_edge_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]
	ids_graph_node_dendrite = \
		list(nx.all_simple_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]

	# Edge lines
	ids_edge_dendrite     = []
	lengths_edges_dendrite = []
	for id_graph_edge in ids_graph_edge_dendrite:
		id    = graph.edges[id_graph_edge]['id']
		loc0  = graph.nodes[id_graph_edge[0]]['loc']
		loc1  = graph.nodes[id_graph_edge[1]]['loc']
		ids_edge_dendrite.append( id )
		lengths_edges_dendrite.append( np.linalg.norm(loc0-loc1) )

	# Tangents
	locations = []
	for id_graph_node in ids_graph_node_dendrite:
		loc = graph.nodes[id_graph_node]['loc']
		locations.append( loc )
	locations = np.array(locations)
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

	return ids_graph_node_dendrite, ids_graph_edge_dendrite, ids_edge_dendrite, lengths_edges_dendrite, tangents, locations



if __name__ == '__main__':
	filename_hdf5 = 'data2/0000000001.hdf5'
	a = CreateGraph(filename_hdf5)
	nx.draw(a.graph, with_labels=True)
	plt.show()

