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
import subprocess as s

import networkx as nx
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import QMainWindow, QTabWidget, QApplication, \
    qApp, QWidget, QHBoxLayout, QVBoxLayout, QLabel, \
    QPushButton, QGraphicsScene, QGraphicsView, QFrame
from PyQt5.QtGui import QFont

from scipy.interpolate import interp1d

import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class Simulation(QMainWindow, PlotCompartmentModelBackend):
	def __init__(self):
		super(Simulation, self).__init__()
		self.initUI2()
		
		# Load data
		data_dir = r'data'
		filename = os.path.join(data_dir, 'compartment_model.pickle')
		with open(filename, 'rb') as f:
			graph    = pickle.load(f)
			self.vertices = pickle.load(f)
			self.head_fs  = pickle.load(f)
			self.neck_fs  = pickle.load(f)

		# Create a compartment model
		self.m = SimulateMembranePotential()
		self.m.tend = 1.5
		self.m.load_spatial_model_graph( graph )
		p1 = self.m.insert_nmdar(123)
		p1['event_time'] = [0.1,0.4,0.6,1.0]
		p2 = self.m.insert_nmdar(101)
		p2['event_time'] = [1.0,1.2,1.4]
		p3 = self.m.insert_nmdar(115)
		p3['event_time'] = [0.9,1.1,1.3,1.45]
		# p = m.insert_current(0)
		
		# Create a new graph.
		num_compartment = graph.number_of_edges()
		c = plt.get_cmap('hsv', 150 )
		self.mesh_actors = []
		for id, edge in enumerate( self.m.graph.edges.values() ):
			if 'faces' in edge:
				faces = edge['faces']
			elif 'id_head' in edge:
				faces = self.head_fs[edge['id_head']]
			elif 'id_neck' in edge:
				faces = self.neck_fs[edge['id_neck']]

			actor = self.plot_mesh(self.vertices, faces, color=c(id)[:3] )
			# actor = self.plot_mesh(self.vertices, faces, color=c(edge['id_compartment'])[:3] )
			bary = obtain_mesh_barycenter(self.vertices, faces)
			#self.plot_text_Billboard(text = str( id ), pos = bary, size = 0.003) # edge['id_compartment']
			#self.removable_actors.AddItem(actor)
			self.renderer.AddActor(actor)
			self.mesh_actors.append( actor )

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
		# self.style.set_node_src_id_dst_id( [self.nodes[0], self.nodes[-1]] )
		self.iren.SetInteractorStyle(self.style)

		# Removable actors
		self.removable_actors = vtk.vtkActorCollection()

		# Update dendritc shaft
		btn = QPushButton(text="Run simulation")
		btn.clicked.connect(self.run_simulation)
		btn.setFont(QFont("Ricty Diminished", 14))
		layout.addWidget(btn)

		self.show()
		#self.iren.Initialize()

	def run_simulation(self):
		print('Run simulation')
		
		#self.m.tend = 2.0
		self.m.run()
		func = interp1d(self.m.t, self.m.y)
		# print('func(1.0)', func(1.0))
		

		im = vtk.vtkWindowToImageFilter()
		im.SetInput(self.vtkWidget.GetRenderWindow())
		im.Update()
		writer = vtk.vtkPNGWriter()
		writer.SetFileName("file.png")
		writer.SetInputData(im.GetOutput())
		writer.Write()


		
		#c = plt.get_cmap('hot', num_compartment )
		c = plt.get_cmap('hot', 100)
		actor_txt_time = []
		for file in glob.glob('imgs/img_*.png'):
  			os.remove(file)
		
		output_dir = 'imgs'
		for id, t in enumerate( np.arange(0.0, self.m.tend, 0.005) ):
			
			y = func(t)
			# print('y ', y)
			for i, actor in enumerate( self.mesh_actors ):
				tmp = int( (y[i] - self.m.p['El']) / 0.05 * 100 )
				if (tmp >= 99):
					tmp = 99
				if (tmp < 0):
					tmp = 0
				actor.GetProperty().SetColor( *c(tmp)[:3] )
			if actor_txt_time != []:
				self.renderer.RemoveActor(actor_txt_time)
			actor_txt_time = self.plot_text_2d(text = "{:.3f} s".format(t))
			self.renderer.AddActor(actor_txt_time)
			self.iren.Initialize()
			
			
			im = vtk.vtkWindowToImageFilter()
			im.SetInput(self.vtkWidget.GetRenderWindow())
			im.Update()
			writer = vtk.vtkPNGWriter()
			fname = os.path.join(output_dir, str(id).zfill(4)+ '.png')
			writer.SetFileName( fname )
			writer.SetInputData(im.GetOutput())
			writer.Write()
			
		
		print('\nConnect images into a video')
		ffname = os.path.join(output_dir, '%04d.png')
		targ_name = 'test_membrane_pot'
		com = ['ffmpeg','-r', '10', '-i', ffname,'-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2','-pix_fmt', 'yuv420p', targ_name+'.mp4']
		print(' '.join(com))
		s.call(com)

		
		'''
		for i in range(self.m.n_compartments):
			plt.plot(self.m.t, self.m.y[i,:],'-', label=str(i))

		plt.legend()
		plt.show()
		'''

if __name__ == '__main__':
	app = QApplication([])
	window = Simulation()
	window.show()
	app.exec_()
	#nx.draw(window.new_graph, with_labels=True)
	#plt.show()


	"""
	self.plot_nodes()
	self.plot_edges()
	self.plot_spines()
	for spine in self.spines:
		for neck, heads in zip( spine['neck'], spine['head'] ):
			if isinstance(heads, list):
				id_neck   = neck['id']
				#print('id_neck ', id_neck)
				container_fs = neck['faces_for_branch']
				
				'''
				for e in self.graph_org.edges( neck['node'] ):
					g_e = self.graph_org.edges[e]
					#print("g_e['path']", g_e['path'])
					paths     = g_e['path'][::20,:]
					radiuses  = g_e['radiuses'][::20]
					normals   = g_e['tangents'][::20,:]
					self.plot_disks( paths, normals, radiuses )

				'''
				
				cols = plt.get_cmap('tab20c', len(container_fs) )# 
				for i in range( len(container_fs) ):
					actor = self.plot_mesh(self.vertices, container_fs[i], color = cols(i)[:3] )
					self.renderer.AddActor(actor)
	
	self.plot_faces_dendrite(self.vertices, self.container_faces)
	self.plot_edges_dendrite(self.locations)

	"""
