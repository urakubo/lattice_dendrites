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

from CreateGraph import CreateGraph
from Utils import *
from DendriticCompartments  import DendriticCompartments
from SpineCompartments      import SpineCompartments
from PlotCompartmentModelBackend import PlotCompartmentModelBackend, Interactor


class Compartmentalization(QMainWindow, PlotCompartmentModelBackend):

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
		self.initial_node_src_id_dst_id = c.find_initial_nodes_dendrite()

		# Dendritic face
		self.vertices, self.faces   = load_stl(fname_mesh)
		#'''
		dend_f        = load_paint(fnames_paint_dendrite[0], self.faces)
		self.dendrite = DendriticCompartments( self.graph, self.vertices, dend_f )
		#'''
		
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
		super(Compartmentalization, self).__init__()
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

		# Show the ids of heads and necks.
		self.delete_billboard()
		for i in range(self.barycenters_head[:,0].shape[0]):
			self.plot_text_Billboard(text = str(i), pos = self.barycenters_head[i,:], size = 0.003)
		for i in range(self.barycenters_neck[:,0].shape[0]):
			self.plot_text_Billboard(text = str(i), color = (0.0,0.0,1.0), pos = self.barycenters_neck[i,:], size = 0.003)

		with open('data_shared.pickle', 'wb') as f:
			pickle.dump(self.graph   , f)
			pickle.dump(self.vertices, f)

		#'''
		print('spine compartments' )
		s  =  SpineCompartments( self.graph, self.style.node_src_id_dst_id, self.vertices, self.head_fs, self.neck_fs, self.barycenters_head )
		s.create()
		pprint.pprint(s.spines)

		with open('data_spines.pickle', 'wb') as f:
			pickle.dump(s.spines     , f)
			pickle.dump(self.head_fs , f)
			pickle.dump(self.neck_fs , f)
		#'''


		#'''
		print('dendrtic compartments' )
		self.dendrite.update_nodes( self.style.node_src_id_dst_id )
		nodes     = self.dendrite.nodes
		locations = [self.graph.nodes[node]['loc'] for node in nodes]
		self.dendrite.create()
		edges     = obtain_path_edges( self.graph, self.style.node_src_id_dst_id )
		
		dendritic_nodes = [{'node': node, 'loc': loc} for node, loc in zip(nodes, locations)]
		dendritic_edges = \
			[{'edge': e, 'faces': f, 'volume': v, 'area': a, 'length': d} \
			for e, f, v, a, d in zip(edges, \
				self.dendrite.split_faces, \
				self.dendrite.split_volumes, \
				self.dendrite.split_areas, \
				self.dendrite.distances)]

		self.plot_faces_dendrite(self.vertices, self.dendrite.split_faces)

		pprint.pprint(dendritic_edges)
		pprint.pprint(dendritic_nodes)
		with open('data_dendrite.pickle', 'wb') as f:
			pickle.dump(dendritic_nodes, f)
			pickle.dump(dendritic_edges, f)
		#'''

		#self.plot_faces_dendrite(vertices, faces)
		self.iren.Initialize()


	def update_dendritc_shaft(self):
		print('Update dendritc shaft')
		self.delete_plots()
		self.delete_billboard()
		for i in range(self.barycenters_head[:,0].shape[0]):
			self.plot_text_Billboard(text = str(i), pos = self.barycenters_head[i,:], size = 0.003)
		for i in range(self.barycenters_neck[:,0].shape[0]):
			self.plot_text_Billboard(text = str(i), color = (0.0,0.0,1.0), pos = self.barycenters_neck[i,:], size = 0.003)

		nodes     = obtain_path_nodes( self.graph,  self.style.node_src_id_dst_id )
		locations = obtain_nodes_location(self.graph, nodes)
		self.plot_edges_dendrite(locations)

		#'''
		self.dendrite.update_nodes( self.style.node_src_id_dst_id )
		self.dendrite.create()
		self.plot_faces_dendrite(self.vertices, self.dendrite.split_faces)
		#'''

		self.iren.Initialize()


if __name__ == '__main__':
    app = QApplication([])
    window = Compartmentalization()
    window.show()
    app.exec_()

