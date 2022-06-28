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
import matplotlib.pyplot as plt


from PyQt5.QtWidgets import QMainWindow, QTabWidget, QApplication, \
    qApp, QWidget, QHBoxLayout, QVBoxLayout, QLabel, \
    QPushButton, QGraphicsScene, QGraphicsView, QFrame
from PyQt5.QtGui import QFont

import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from CreateCompartmentModel import CreateGraph, obtain_edges_dendrite, load_stl, load_paint, fill_hole
from PlotCompartmentModelBackend import PlotCompartmentModelBackend, Interactor


class PlotCompartmentModel(QMainWindow, PlotCompartmentModelBackend):

	def _loader(self):
		annot_dir = r"C:\Users\uraku\Desktop\LatticeMicrobes\220610ReconstMorph\annots\dend4_220610_without_spine"
		annot_dir = annot_dir.replace('/', os.sep)
		fname_mesh   = os.path.join(annot_dir, "surfaces","whole","0000000001.stl")
		fname_hdf5   = os.path.join(annot_dir, "skeletons","whole","0000000001.hdf5")
		fname_paint  = os.path.join(annot_dir, "paint","0000000001-*.pickle")
		fnames_paint = glob.glob(fname_paint)

		### Skeleton
		c          = CreateGraph(fname_hdf5)
		self.graph = c.graph
		self.node_src_id_dst_id = c.node_src_id_dst_id

		### Dendritic face
		vertices, faces   = load_stl(fname_mesh)
		'''
		self.v_org , self.f_org , self.fcenter_org , self.farea_org = load_paint(fnames_paint[0], vertices, faces)
		self.v_fill, self.f_fill, self.fcenter_fill, self.farea_fill, volume = fill_hole(self.v_org, self.f_org)
		'''
		
		### Spine-head face
		annot_dir     = r"C:\Users\uraku\Desktop\LatticeMicrobes\220610ReconstMorph\annots\dend4_220610_head"
		annot_dir     = annot_dir.replace('/', os.sep)
		fname_mesh   = os.path.join(annot_dir, "surfaces","whole","0000000001.stl")
		fname_paint   = os.path.join(annot_dir, "paint","0000000001-*.pickle")
		fnames_paint2 = glob.glob(fname_paint)
		#vertices2, faces2   = load_stl(fname_mesh)
		self.v_spine = []
		self.f_spine = []
		for fname in fnames_paint2:
			#print('fname ', fname)
			v_spine, f_spine, _, _ = load_paint(fname, vertices, faces)
			self.v_spine.append(v_spine)
			self.f_spine.append(f_spine)

		### How do I obtain spine-neck faces?



	def __init__(self):
		super(PlotCompartmentModel, self).__init__()
		self._loader()
		self.initUI()

		self.plot_nodes()
		self.plot_edges()
		self.refresh_plots()

	def refresh_plots(self):
		print('Refresh plots')
		self.delete_plots()

		# Obtain dendrite and plot dendritic skeleton (blue line)
		ids_graph_node_dendrite, ids_graph_edge_dendrite, ids_edge_dendrite, \
				lengths_edge_dendrite, tangents, locations \
				= obtain_edges_dendrite( self.graph, self.style.node_src_id_dst_id )

		self.plot_edges_dendrite(locations)
		
		# Plot dendrite
		# volumes, areas = self.plot_faces_dendrite(tangents, locations)


		# Plot spines
		c = plt.get_cmap('cool', len(self.v_spine)+1 ) # hsv, flag
		for i, (v, f) in enumerate( zip(self.v_spine, self.f_spine) ):
			#print("f ", f)
			actor = self.plot_mesh(v, f, color = c(i)[:3] )
			self.renderer.AddActor(actor)


		self.iren.Initialize()

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


	def initUI(self):
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
		self.style.set_node_src_id_dst_id( self.node_src_id_dst_id )
		self.iren.SetInteractorStyle(self.style)
		# self.style.node_src_id_dst_id = self.set_node_src_id_dst_id

		# Removable actors
		self.removable_actors = vtk.vtkActorCollection()

		#self.viewer.setCameraPosition(distance=10)
		btn = QPushButton(text="Refresh screen")
		btn.clicked.connect(self.refresh_plots)
		btn.setFont(QFont("Ricty Diminished", 14))
		layout.addWidget(btn)

		self.show()
		#self.iren.Initialize()


if __name__ == '__main__':
    app = QApplication([])
    window = PlotCompartmentModel()
    window.show()
    app.exec_()

