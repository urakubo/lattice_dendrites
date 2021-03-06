import numpy as np
import networkx as nx
import trimesh
import pymeshfix
import itertools


class CloseMesh():
	"""Create the closed mesh from a open mech using the pymeshfix. Utility functions are also implemented.

	Args:
		vertices (numpy[float]): Vertices of a triangular mesh. vertices.shape is (n, 3).
		faces (numpy[int]): Vertical array of three vertices that specify a triangle. faces.shape is (m, 3).

	Returns:
		(pyLD.CloseMesh): CloseMesh object that has the follwing instances:
		- unclosed_area (float): Summed area of input faces.
		- vertices (numpy[float]): Input vertices
		- faces (numpy[float]): Input faces
		- mesh (trimesh object): Closed mesh
		- volume (float): Volume of the closed mesh.
	"""
	def __init__(self, vertices, faces):

		self.mesh0 = trimesh.Trimesh(vertices, faces)
		self.unclosed_area = self.mesh0.area

		self.vertices = vertices
		self.faces    = faces

		part_mesh = pymeshfix.MeshFix(vertices, faces)
		part_mesh.repair()
		self.mesh = trimesh.Trimesh(part_mesh.v, part_mesh.f)
		self.mesh.merge_vertices()
		self.mesh.remove_degenerate_faces()
		self.mesh.remove_duplicate_faces()
		self.volume    = self.mesh.volume

	def obtain_vertices(self):
		"""Obtain vertices of the closed mesh.
		Returns:
			- vertices (numpy[float]): vertices.shape is (n, 3). 
		"""
		return np.array(self.mesh.vertices)

	def obtain_faces(self):
		"""Obtain faces of the closed mesh.
		Returns:
			- vertices (numpy[int]): vertices.shape is (n, 3). 
		"""
		return np.array(self.mesh.faces)

	def obtain_fcenters(self):
		"""Obtain the center coordinates of faces in the closed mesh.
		Returns:
			- vertices (numpy[int]): vertices.shape is (n, 3). 
		"""
		return np.array(self.mesh.triangles_center)

	def obtain_crossing_edges(self, graph):
		"""Obtain the edges that cross the closed volume if the networkx graph has 'loc' in the nodes.
		Returns:
			- crossing_edges (list[(int, int)]): List of the edges that cross the closed volume.
		"""
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
		"""Obtain the length of the edge of a networkx graph.

		Args:
			graph (Networkx graph object) : It must contain node['loc'] and edge['path'], both of which specify the location.

		Returns:
			- (numpy.float): the total length of edges.
		"""
		return sum( [self.obtain_inside_length_edge( edge ) for edge in graph.edges.values()] )

	def obtain_lengths_inside(self, graphs):
		"""Obtain the length of the edge of networkx graphs.

		Args:
			graphs (list[Networkx graph object]) : They must contain node['loc'] and edge['path'], both of which specify the location.

		Returns:
			- (list[numpy.float]): the total lengths of graphs.
		"""
		return [self.obtain_length_inside( graph ) for graph in graphs ]

	def assign_mesh_for_each_graph(self, sub_gs):
		paths_enclosed = [self._pick_up_5_enclosed_points(sub_g) for sub_g in sub_gs]
		# print('paths_enclosed ', paths_enclosed)
		nums_points    = [path.shape[0] for path in paths_enclosed]
		num_branches   = len(nums_points)
		fcenters = self.mesh0.triangles_center
		dists = []
		for num_points, path in zip(nums_points, paths_enclosed):
			dist = np.array( [np.linalg.norm((fcenters - path[i,:]), axis=1) for i in range(num_points)] )
			#print('Before aaggregation dist.shape', dist.shape)
			dist = np.min( dist, axis = 0)
			#print('After aaggregation  dist.shape', dist.shape)
			dists.append(dist)
		dists = np.array(dists)
		ids_face_for_branch = np.argmin(dists, axis=0)

		faces_for_branch  = [ self.faces[ids_face_for_branch==i,:] for i in range(num_branches) ]
		meshes_for_branch = [ trimesh.Trimesh(self.vertices, f)    for f in faces_for_branch    ]
		fareas_for_branch = [ m.area                               for m in meshes_for_branch   ]
		# print('closed_ids  ', closed_ids)
		# print('dists.shape ', dists.shape)
		# print('closed_ids.shape', closed_ids.shape)
		return faces_for_branch, fareas_for_branch


	def obtain_open_mesh_pertial_area(self, ids):
		vertices = self.mesh0.vertices
		faces    = self.mesh0.faces[ids]
		mesh_partial = trimesh.Trimesh(vertices, faces)
		return mesh_partial.area


	def _pick_up_5_enclosed_points(self, sub_g):
		p_enclosed = np.empty([0,3])
		for edge in sub_g.edges.values():
			p = edge['path']
			p = p[self.mesh.contains(p) == True,:]
			p_enclosed = np.vstack([p_enclosed, p])
		print('Num of path points in the target volume (min > 1)', p_enclosed.shape[0])
		if p_enclosed.shape[0] < 2:
			raise TypeError('Error! Too small number of path points are enclosed.')
		elif p_enclosed.shape[0] > 5:
			splitted   = np.array_split(p_enclosed, 5, 0)
			p_enclosed = np.array( [s[0,:] for s in splitted] )
		return p_enclosed

	def obtain_nodes_inside(self, graph):
		nodes_inside = []
		for id_node in graph.nodes.keys():
			loc0 = graph.nodes[id_node]['loc'].reshape([1,-1])
			if self.mesh.contains(loc0):
				nodes_inside.append(id_node)
		return nodes_inside


	def obtain_inside_length_edge(self, edge):
		"""Obtain the length of the edge in the closed volume.

		Args:
			edge (the edge of a networkx graph: Vertices of a trace. edge['path'].shape is (n, 3).

		Returns:
			- ength_sub_path (numpy.float): the length of the edge.
		"""
		path    = edge['path']
		# points ((n, 3) float)
		# return (n, ) bool
		enclosed = self.mesh.contains(path)
		sub_path = path[enclosed == True,:]
		
		diff_sub_path   = np.diff(sub_path, axis=0)
		length_sub_path = np.sum(np.linalg.norm(diff_sub_path, axis=1))
		#print('length_sub_path ', length_sub_path)
		return length_sub_path

