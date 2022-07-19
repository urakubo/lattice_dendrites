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

from vtk_version.Utils import *
from vtk_version.PlotCompartmentModelBackend import PlotCompartmentModelBackend, Interactor


class PreprocessSimulation(QMainWindow, PlotCompartmentModelBackend):

	def _loader(self):
		data_dir = r'C:\Users\uraku\Desktop\LatticeMicrobes\sim_membranepot\vtk_version'
		
		filename = os.path.join(data_dir, 'data_shared.pickle')
		with open(filename, 'rb') as f:
			self.graph = pickle.load(f)
			self.vertices = pickle.load(f)

		filename = os.path.join(data_dir, 'data_spines.pickle')
		with open(filename, 'rb') as f:
			self.spines   = pickle.load(f)
			self.head_fs  = pickle.load(f)
			self.neck_fs  = pickle.load(f)

		filename = os.path.join(data_dir, 'data_dendrite.pickle')
		with open(filename, 'rb') as f:
			self.dendritic_nodes = pickle.load(f)
			self.dendritic_edges = pickle.load(f)

		self.nodes           = [n['node'] for n in self.dendritic_nodes]
		self.locations       = obtain_nodes_location(self.graph, self.nodes)
		self.container_faces = [n['faces'] for n in self.dendritic_edges]

	def __init__(self):
		super(PreprocessSimulation, self).__init__()
		self._loader()
		self.initUI2()

		self.plot_nodes()
		self.plot_edges()
		#self.plot_spines()
		
		
		for spine in self.spines:
			for neck, heads in zip( spine['neck'], spine['head'] ):
				if isinstance(heads, list):
					id_neck   = neck['id']
					#print('id_neck ', id_neck)
					container_fs = neck['faces_for_branch']
					cols = plt.get_cmap('hsv', len(container_fs) )
					
					for e in self.graph.edges( neck['node'] ):
						g_e = self.graph.edges[e]
						#print("g_e['path']", g_e['path'])
						paths     = g_e['path'][::20,:]
						radiuses  = g_e['radiuses'][::20]
						normals   = g_e['tangents'][::20,:]
						self.plot_disks( paths, normals, radiuses )
					
					for i in range( len(container_fs) ):
						actor = self.plot_mesh(self.vertices, container_fs[i], color = cols(i)[:3] )
						self.renderer.AddActor(actor)
		
		
		self.plot_faces_dendrite(self.vertices, self.container_faces)
		self.plot_edges_dendrite(self.locations)
		self.iren.Initialize()


	def initUI2(self):
		self.setGeometry(0, 0, 700, 900) 
		centerWidget = QWidget()
		self.setCentralWidget(centerWidget)

		layout = QVBoxLayout()
		centerWidget.setLayout(layout)

		self.frame = QFrame()
		self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
		layout.addWidget(self.vtkWidget, 1)

		# Vtk root
		self.renderer = vtk.vtkRenderer()
		self.renderer.SetBackground(0.1, 0.2, 0.4)
		self.vtkWidget.GetRenderWindow().AddRenderer(self.renderer)
		self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

		# style = vtk.vtkInteractorStyleTrackballCamera()
		self.style = Interactor()
		self.style.set_renderer( self.renderer )
		self.style.set_node_src_id_dst_id( [self.nodes[0], self.nodes[-1]] )
		self.iren.SetInteractorStyle(self.style)

		# Removable actors
		self.removable_actors = vtk.vtkActorCollection()

		self.show()
		#self.iren.Initialize()



if __name__ == '__main__':
    app = QApplication([])
    window = PreprocessSimulation()
    window.show()
    app.exec_()

