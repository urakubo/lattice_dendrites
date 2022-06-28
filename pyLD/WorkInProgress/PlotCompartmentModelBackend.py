#
#


import numpy as np
import matplotlib.pyplot as plt

import sys
import vtk

#import networkx as nx
#import itertools
#from stl import mesh
from CreateCompartmentModel import CreateGraph, obtain_edges_dendrite, load_stl, load_paint, fill_hole


class Interactor(vtk.vtkInteractorStyleTrackballCamera):
	def __init__(self, parent=None):
		self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
		self.renderer           = None
		self.node_src_id_dst_id = None
		self.pickuppable_actors = vtk.vtkActorCollection()

	def set_renderer(self, renderer):
		self.renderer = renderer

	def set_node_src_id_dst_id(self, node_src_id_dst_id):
		self.node_src_id_dst_id = node_src_id_dst_id
		print('node_src_id_dst_id: ', self.node_src_id_dst_id)

	def leftButtonPressEvent(self, obj, event):

		clickPos = self.GetInteractor().GetEventPosition()
		#print(f'left button down: {clickPos}')

		picker = vtk.vtkPropPicker()
		# picker.Pick(clickPos[0], clickPos[1], 0, self.renderer)
		picker.PickProp(clickPos[0], clickPos[1], self.renderer, self.pickuppable_actors)
		info = picker.GetActor()
		# Create a new actor
		# print('info: ', info)
		if (info is not None) and hasattr(info, 'key'):
			summary = info.GetProperty().GetInformation()
			name    = info.key
			#print('summary ', summary)
			#print('picker.GetActor().key   ', name )
			if 'node' in name:
				id = int(name.split('_')[1])
				self.node_src_id_dst_id = [self.node_src_id_dst_id[1], id]
				print(self.node_src_id_dst_id)


		self.OnLeftButtonDown()


class PlotCompartmentModelBackend():

	def delete_plots(self):
		for i in range( self.removable_actors.GetNumberOfItems() ):
			select_actor = self.removable_actors.GetItemAsObject(i)
			self.renderer.RemoveActor(select_actor)
		self.removable_actors.RemoveAllItems()

	def plot_mesh(self, vertices, faces, color = (1.0, 1.0, 1.0) ):
	
		if faces is None:
			actor = vtk.vtkActor()
			return actor
	
		polydata = self.CreatePolyData(vertices, faces)
		# mapper
		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputData( polydata )
		# actor
		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		actor.GetProperty().SetColor(*color)
		actor.GetProperty().SetOpacity(0.5)
		return actor


	def CreatePolyData( self, pts, faces ):
		"""
		Creates vtkPolyData from vertices and faces

		pts numpy.array: Nx3 array of vertices
		faces numpy.array: Mx3 array of faces

		Return vtkPolyData
		"""
		(nv,mv) = pts.shape
		(nf,mf) = faces.shape
		cells = vtk.vtkCellArray()
		for j in range(nf):
			cell = vtk.vtkTriangle()
			cell.GetPointIds().SetNumberOfIds(3)
			cell.GetPointIds().SetId( 0, faces[j,0] )
			cell.GetPointIds().SetId( 1, faces[j,1] )
			cell.GetPointIds().SetId( 2, faces[j,2] )
			cells.InsertNextCell( cell )


		points = vtk.vtkPoints()
		points.SetNumberOfPoints(nv)
		for j in range(nv):
			points.SetPoint( j, pts[j,0], pts[j,1], pts[j,2] )

		new_mesh = vtk.vtkPolyData()
		new_mesh.SetPoints( points )
		new_mesh.SetPolys( cells )
		new_mesh.BuildCells()	

		return new_mesh


	def plot_nodes(self):
		for id, node in self.graph.nodes(data=True):
			pos    = node['loc']
			# Generate a sphere for the node
			sphere = vtk.vtkSphereSource()
			sphere.SetCenter(*pos)
			sphere.SetRadius(0.05)

			# Set up a mapper for the node
			mapper = vtk.vtkPolyDataMapper()
			mapper.SetInputConnection(sphere.GetOutputPort())

			# Set up an actor for the node
			actor = vtk.vtkActor()
			actor.key = "node_"+str(id)
			

			actor.GetProperty().SetColor(1.0,0.0,0.0) # Yellow
			actor.SetMapper(mapper)

			self.style.pickuppable_actors.AddItem(actor)
			self.renderer.AddActor(actor)


	def _plot_line(self, pos, color = (1.0,1.0,1.0), width=1.0):
		pnum = pos.shape[0]
		# print('pnum: ', pnum)
		points = vtk.vtkPoints()
		for i in range(pnum):
			# print('pos[i,:] ', pos[i,:])
			points.InsertNextPoint(pos[i,:])
		polyLine = vtk.vtkPolyLine()
		polyLine.GetPointIds().SetNumberOfIds(pnum)
		for i in range(pnum):
			polyLine.GetPointIds().SetId(i, i)
		cells = vtk.vtkCellArray()
		cells.InsertNextCell(polyLine)
		polyData = vtk.vtkPolyData()
		polyData.SetPoints(points)
		polyData.SetLines(cells)
		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputData(polyData)
		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		actor.GetProperty().SetColor(*color)
		actor.GetProperty().SetLineWidth(width)
		return actor

	def plot_edges(self):
		for u,v, data in self.graph.edges(data=True):
			# Draw lines
			pos   = data['path']
			id    = data['id']
			actor = self._plot_line(pos)
			actor.key = "edge_"+str(id)
			self.renderer.AddActor(actor)

	def plot_edges_dendrite(self, locations):
		# Edge line
		lengths = []
		color_blue  =(0.0, 0.0, 1.0)
		actor = self._plot_line(locations, color=color_blue, width=5.0)
		actor.key = "edge_dend_"+str(id)
		self.removable_actors.AddItem(actor)
		self.renderer.AddActor(actor)

	def plot_grid(self):
		axes = vtk.vtkAxesActor()
		axes_widget = vtk.vtkOrientationMarkerWidget()
		axes_widget.SetOutlineColor(1, 1, 1)
		axes_widget.SetOrientationMarker(axes)
		axes_widget.SetInteractor(inter)
		axes_widget.SetViewport(0., 0., 0.4, 0.4)
		axes_widget.EnabledOn()
		axes_widget.InteractiveOff()

	def plot_faces_dendrite(self, tangents, locations_):

		locations  = np.delete( locations_, [0, -1], 0)

		# Obtain faces
		ids_f_org_part  = {}
		ids_f_fill_part = {}
		for i in range(len(tangents)):
			tangent  = tangents[i]
			location = locations[i]
			ids_f_fill_part[i] = ( np.dot(self.fcenter_fill-location, tangent) > 0 )
			ids_f_org_part[i]  = ( np.dot(self.fcenter_org -location, tangent) > 0 )


		# Surface areas, volumes and visualzation
		c = plt.get_cmap('hsv', len(tangents)+1 ) # hsv, flag
		volumes  = []
		areas    = []
		for i in range(len(tangents)+1):
			if i == 0:
				i_part_f_fill = ids_f_fill_part[i]
				i_part_f_org  = ids_f_org_part[i]
			elif i < len(tangents):
				i_part_f_fill = ids_f_fill_part[i] * np.logical_not(ids_f_fill_part[i-1])
				i_part_f_org  = ids_f_org_part[i]  * np.logical_not(ids_f_org_part[i-1])
			elif i == len(tangents):
				i_part_f_fill = np.logical_not(ids_f_fill_part[i-1])
				i_part_f_org  = np.logical_not(ids_f_org_part[i-1])
			else:
				print('Error.')
				break
			part_f_fill = self.f_fill[i_part_f_fill]
			area        = np.sum( self.farea_org[i_part_f_org] )
			vertices, faces, fcenter, farea, volume = fill_hole(self.v_fill, part_f_fill)
			actor = self.plot_mesh(vertices, faces, color=c(i)[:3] )
			self.removable_actors.AddItem(actor)
			self.renderer.AddActor(actor)
			volumes.append(volume)
			areas.append(area)

		return volumes, areas
		#


