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
		self.initial_node_src_id_dst_id = c.node_src_id_dst_id

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
		
		# Reorder spine heads and necks
		self.barycenters_head   = np.array( [obtain_mesh_barycenter(self.vertices, f) for f in self.head_fs] )
		self.barycenters_neck   = np.array( [obtain_mesh_barycenter(self.vertices, f) for f in self.neck_fs] )
		
		order_head = np.argsort( self.barycenters_head[:,2] )
		self.barycenters_head = self.barycenters_head[order_head,:]
		self.head_fs          = [ self.head_fs[id] for id in order_head.tolist() ]
		order_neck = np.argsort( self.barycenters_neck[:,2] )
		self.barycenters_neck = self.barycenters_neck[order_neck,:]
		self.neck_fs          = [ self.neck_fs[id] for id in order_neck.tolist() ]


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

		s  =  SpineCompartments( self.graph, self.style.node_src_id_dst_id, self.vertices, self.head_fs, self.neck_fs, self.barycenters_head )
		s.create()

		# Show the ids of heads and necks.
		self.delete_billboard()
		for i in range(self.barycenters_head[:,0].shape[0]):
			self.plot_text_Billboard(text = str(i), pos = self.barycenters_head[i,:], size = 0.003)
		for i in range(self.barycenters_neck[:,0].shape[0]):
			self.plot_text_Billboard(text = str(i), color = (0.0,0.0,1.0), pos = self.barycenters_neck[i,:], size = 0.003)


		print('graphs_spine ' )
		pprint.pprint(s.spines)


		with open('data_spines.pickle', 'wb') as f:
			pickle.dump(self.graph   , f)
			pickle.dump(s.spines     , f)
			pickle.dump(self.vertices, f)
			pickle.dump(self.head_fs , f)
			pickle.dump(self.neck_fs , f)

		# ids     = obtain_path_nodes( self.graph, self.style.node_src_id_dst_id )
		# lengths = obtain_distance_between_edges_dendrite( self.graph, self.style.node_src_id_dst_id )
		# lengths_accumulated = [0.0] + list(itertools.accumulate(lengths))


		'''
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
		
		# spines = obtain_graphs_spines_along_dendrite( self.graph, self.style.node_src_id_dst_id )
		# print('graphs_spine ', graphs_spine)

		#for i, i_node in enumerate( ids_node_terminal ):
		#	self.plot_text_Billboard(text = str(i_node), pos =  self.graph.nodes[i_node]['loc'])

		self.iren.Initialize()


if __name__ == '__main__':
    app = QApplication([])
    window = PlotCompartmentModel()
    window.show()
    app.exec_()

