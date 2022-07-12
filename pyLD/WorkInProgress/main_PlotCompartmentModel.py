#
#
#
import os, sys, glob
import numpy as np
import itertools
from pathlib import Path
import h5py
import pprint
import pickle
import trimesh

import networkx as nx
import matplotlib.pyplot as plt


from PyQt5.QtWidgets import QMainWindow, QTabWidget, QApplication, \
    qApp, QWidget, QHBoxLayout, QVBoxLayout, QLabel, \
    QPushButton, QGraphicsScene, QGraphicsView, QFrame
from PyQt5.QtGui import QFont

import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from CreateCompartmentModel import *
from PlotCompartmentModelBackend import PlotCompartmentModelBackend, Interactor


class PlotCompartmentModel(QMainWindow, PlotCompartmentModelBackend):

	def _loader(self):

		# 
		annot_dir = r"C:\Users\uraku\Desktop\LatticeMicrobes\220610ReconstMorph\annots\dend4_220610_without_spine"
		annot_dir = annot_dir.replace('/', os.sep)
		fname_mesh   = os.path.join(annot_dir, "surfaces","whole","0000000001.stl")
		fname_hdf5   = os.path.join(annot_dir, "skeletons","whole","0000000001.hdf5")
		fnames_paint_dendrite = os.path.join(annot_dir, "paint","0000000001-*.pickle")
		fnames_paint_dendrite = glob.glob(fnames_paint_dendrite)

		# Skeleton
		c          = CreateGraph(fname_hdf5)
		self.graph = c.graph
		self.node_src_id_dst_id = c.node_src_id_dst_id

		# Dendritic face
		self.vertices, self.faces   = load_stl(fname_mesh)
		'''
		self.v_org , self.f_org , self.fcenter_org , self.farea_org = load_paint(fnames_paint_dendrite[0], self.vertices, self.faces)
		self.v_fill, self.f_fill, self.fcenter_fill, self.farea_fill, volume = fill_hole(self.v_org, self.f_org)
		'''
		
		# Spine-head face
		annot_dir     = r"C:\Users\uraku\Desktop\LatticeMicrobes\220610ReconstMorph\annots\dend4_220610_head"
		annot_dir     = annot_dir.replace('/', os.sep)
		fname_mesh    = os.path.join(annot_dir, "surfaces","whole","0000000001.stl")
		fnames_paint_head = os.path.join(annot_dir, "paint","0000000001-*.pickle")
		fnames_paint_head = glob.glob(fnames_paint_head)
		self.head_fs = load_paints(fnames_paint_head, self.faces)

		# Spine-neck faces
		fnames_paint = fnames_paint_dendrite + fnames_paint_head
		self.neck_fs = load_paint_subtraction(fnames_paint, self.faces)

	def __init__(self):
		super(PlotCompartmentModel, self).__init__()
		self._loader()
		self.initUI()

		self.plot_nodes()
		self.plot_edges()
		self.plot_spines()

		self.dendritic_node_has = []
		self.ids_node_terminal  = []
		self.update_dendritc_shaft()
		

	def create_model(self):
		print('Create model')
		# ids     = obtain_ids_graph_node_dendrite( self.graph, self.style.node_src_id_dst_id )
		# lengths = obtain_distance_between_edges_dendrite( self.graph, self.style.node_src_id_dst_id )
		# lengths_accumulated = [0.0] + list(itertools.accumulate(lengths))

		# Obtain spine-heads associated with nodes

		# Preparation
		spines, ids_node_terminal = obtain_graphs_spines_along_dendrite( self.graph, self.style.node_src_id_dst_id )
		barycenters_head   = np.array( [obtain_mesh_barycenter(self.vertices, f) for f in self.head_fs] )
		barycenters_neck   = np.array( [obtain_mesh_barycenter(self.vertices, f) for f in self.neck_fs] )
		adjacents = AdjacentFaces( self.vertices, self.head_fs, self.neck_fs )

		# Find spine heads with skeleons
		for spine in spines:
			ids_head  = []
			terminals = []
			lengths   = []
			volumes   = []
			areas     = []
			# crossing_edges = []
			for id in spine['terminals']:
				loc_terminal = spine['graph'].nodes[id]['loc']
				dists  = np.linalg.norm(barycenters_head-loc_terminal, axis=1)
				id_head = np.argmin(dists)
				if id_head in ids_head:
					print('Warning: the selected spine is again selected. Ignored.')
				else:
					ids_head.append(id_head)
					terminals.append(id)
					head = ClosedMesh( self.vertices, self.head_fs[id_head] )
					lengths.append( head.obtain_inside_length_graph( spine['graph']) )
					volumes.append( head.volume )
					areas.append( head.unclosed_area )
					# crossing_edges.append( head.obtain_crossing_edges( spine['graph'] ) )
					
			spine['ids_head']     = ids_head
			spine['terminals']    = terminals
			spine['head_lengths'] = lengths
			spine['head_volumes'] = volumes
			spine['head_areas']   = areas
			# spine['head_crossing_edges'] = crossing_edges


		# Find spine necks
		for spine in spines:
			ids_neck = []
			for id_head in spine['ids_head']:
				id_neck = adjacents.obtain_necks_connecting_to_head( id_head )
				if len(id_neck) >= 2:
					print('Warning: more than two necks are connected to a spine. Only one is connected.')
					id_neck = id_neck[0]
				elif len(id_neck) == 0:
					print('Warning: no spine neck is connected to a spine.')
					id_neck = None
				ids_neck.append(id_neck) # id_neck[0]
			spine['ids_neck'] = ids_neck

		print('graphs_spine ' )
		pprint.pprint(spines)


		# Calculate paths
		for spine in spines:
			ids_head = spine['ids_head']
			ids_neck = spine['ids_neck']
			
			# if there are no shared branches:
			if len(ids_neck) == len(set(ids_neck)): # It accept ids_neck = [None]
				for id_head, id_neck in zip(ids_head, ids_neck):
					pass
			# if it has shared branches
			else :
				ids, counts = np.unique(ids_neck, return_counts=True)



			for id_head in ids_head:
				print( 'spine head id: ', id_head )
				head = ClosedMesh( self.vertices, self.head_fs[id_head] )
				head_length = [head.obtain_inside_length( edge ) for edge in spine['graph'].edges.values()]
				print('head_length ', head_length )



#		print('graphs_spine ' )
#		pprint.pprint(spines)


		#for i, (id_head, id_neck) in enumerate( zip(ids_head, ids_neck) ):
		#	spine_head = ClosedMesh( self.vertices, self.faces )
		#	contains(points)

		# Obtain spine-necks associated with spine-heads
		#obtain_necks_connecting_to_each_heads()


		'''		
		dendritic_node_has = obtain_ids_of_terminal_nodes_from_dendrite( self.graph, self.style.node_src_id_dst_id )

		print('ids ', ids)
		print('lengths_accumulated ', lengths_accumulated)
		
		ids_locs = list(self.graph.nodes(data='loc'))
		print('ids_locs ',ids_locs)
		pos = { id: (loc[1], loc[2]) for id, loc in ids_locs}
		print(pos)
		nx.draw_networkx(self.graph, pos=pos) # , node_color="c"
		plt.show()

		### Obtain adjacency
		heads_boundary_vertices = [obtain_boundary_vertices(self.vertices, f) for f in self.head_fs]
		necks_boundary_vertices = [obtain_boundary_vertices(self.vertices, f) for f in self.neck_fs]
		heads_bind_to_neck = []
		for i, neck_vertices in enumerate(necks_boundary_vertices):
			adjacents = []
			for j, head_vertices in enumerate(heads_boundary_vertices):
				#print('len(ref_boundary_vertices) ', len(ref_boundary_vertices))
				if neck_vertices & head_vertices != set():
					adjacents.append(j)
			heads_bind_to_neck.append(adjacents)
		print('heads_bind_to_neck ', heads_bind_to_neck )
		'''
		'''
		print('edge ids len : ', len(ids_edge_dendrite))
		print('lengths  len : ', len(lengths_edge_dendrite))
		print('volumes  len : ', len(volumes))
		print('areas    len : ', len(areas))

		with open('data.pickle', 'wb') as f:
			pickle.dump(ids_graph_node_dendrite,f)
			pickle.dump(ids_graph_edge_dendrite,f)
			pickle.dump(ids_edge_dendrite      ,f)
			pickle.dump(lengths_edge_dendrite  ,f)
			pickle.dump(volumes  ,f)
			pickle.dump(areas    ,f)
			pickle.dump(locations,f)
			pickle.dump(tangents ,f)
		'''

	def update_dendritc_shaft(self):
		print('Update dendritc shaft')
		self.delete_plots()

		# self.plot_text_2d()
		# Obtain and plot the skeleton (blue line) of dendrite

		locations, tangents = obtain_locs_tangents_from_dendrite( self.graph, self.style.node_src_id_dst_id )
		self.plot_edges_dendrite(locations)
		#volumes, areas = self.plot_faces_dendrite(tangents, locations)

		# Plot revised node IDs
		self.delete_billboard()
		
		'''
		dendritic_node_has = obtain_ids_of_terminal_nodes_from_dendrite( self.graph, self.style.node_src_id_dst_id )
		ids_node_terminal  = list(itertools.chain.from_iterable( dendritic_node_has.values() ))
		'''
		
		graphs_spine, ids_node_terminal = obtain_graphs_spines_along_dendrite( self.graph, self.style.node_src_id_dst_id )
		#print('graphs_spine ', graphs_spine)

		for i, i_node in enumerate( ids_node_terminal ):
			self.plot_text_Billboard(text = str(i), pos =  self.graph.nodes[i_node]['loc'])


		self.iren.Initialize()


if __name__ == '__main__':
    app = QApplication([])
    window = PlotCompartmentModel()
    window.show()
    app.exec_()

