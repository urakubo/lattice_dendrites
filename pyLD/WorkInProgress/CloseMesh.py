import numpy as np
import networkx as nx
import trimesh
import pymeshfix
import itertools



class CloseMesh():
	def __init__(self, vertices, faces):

		mesh0 = trimesh.Trimesh(vertices, faces)
		self.unclosed_area = mesh0.area

		part_mesh = pymeshfix.MeshFix(vertices, faces)
		part_mesh.repair()
		self.mesh = trimesh.Trimesh(part_mesh.v, part_mesh.f)
		self.mesh.merge_vertices()
		self.mesh.remove_degenerate_faces()
		self.mesh.remove_duplicate_faces()
		self.volume    = self.mesh.volume
	
	def obtain_vertices(self):
		return np.array(self.mesh.vertices)

	def obtain_faces(self):
		return np.array(self.mesh.faces)

	def obtain_fcenters(self):
		return np.array(self.mesh.triangles_center)

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

