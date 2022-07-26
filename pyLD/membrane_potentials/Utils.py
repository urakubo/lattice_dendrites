import numpy as np
import networkx as nx
import pickle
import gzip
import trimesh
import glob
import itertools
import copy


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


def load_paint(filename_paint, faces):
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
		sub_face = load_paint(filename_paint, faces)
		multiple_faces.extend( _separate_unconnected_faces(sub_face) )
	return multiple_faces


def obtain_nodes_location(graph, nodes):
	locations = [graph.nodes[node]['loc'] for node in nodes]
	locations = np.array(locations)
	return locations


def obtain_path_nodes(graph, src_id_dst_id):
	return list(nx.all_simple_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]


def obtain_path_edges(graph, src_id_dst_id):
	return list(nx.all_simple_edge_paths(graph, src_id_dst_id[0], src_id_dst_id[1]))[0]

