#
#
import numpy as np
import matplotlib.pyplot as plt
import sys

from PyQt5.QtWidgets import QMainWindow, QTabWidget, QApplication, \
    qApp, QWidget, QHBoxLayout, QVBoxLayout, QLabel, \
    QPushButton, QGraphicsScene, QGraphicsView, QFrame
from PyQt5.QtGui import QFont
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class Interactor(vtk.vtkInteractorStyleTrackballCamera):
	def __init__(self, parent=None):
		self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
		self.renderer              = None
		self.node_src_id_dst_id    = None
		self.pickuptable_actors    = vtk.vtkActorCollection()
		self.three_d_text_actors   = []
		self.text_Billboard_actors = []


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
		picker.PickProp(clickPos[0], clickPos[1], self.renderer, self.pickuptable_actors)
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

	'''
	def leftButtonReleaseEvent(self, obj, event):
		camera = self.renderer.GetActiveCamera()
		#x, y, z = camera.GetDirectionOfProjection() # Get the vector in the direction from the camera position to the focal point.
		#print('x, y, z ', x, y, z)
		orientation = camera.GetOrientationWXYZ()
		
		# RotateZ = np.rad2deg(np.arctan(y/x) )
		
		
		for act_txt in self.three_d_text_actors:
			# act_txt.SetOrientation(orientation[1] ,orientation[0] ,orientation[2] )
			# specified as X, Y and Z, and performed as RotateZ, RotateX, and finally RotateY.
			act_txt.RotateWXYZ( *orientation )
		self.OnLeftButtonUp()
	'''

class PlotCompartmentModelBackend():


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
		self.style.set_node_src_id_dst_id( self.initial_node_src_id_dst_id )
		self.iren.SetInteractorStyle(self.style)
		# self.style.node_src_id_dst_id = self.set_node_src_id_dst_id

		# Removable actors
		self.removable_actors = vtk.vtkActorCollection()

		# Update dendritc shaft
		btn = QPushButton(text="Update dendritic shaft")
		btn.clicked.connect(self.update_dendritc_shaft)
		btn.setFont(QFont("Ricty Diminished", 14))
		layout.addWidget(btn)

		# Create multicompartment model
		btn = QPushButton(text="Create multicompartment model")
		btn.clicked.connect(self.create_model)
		btn.setFont(QFont("Ricty Diminished", 14))
		layout.addWidget(btn)

		self.show()
		#self.iren.Initialize()

		"""
		camera = self.renderer.GetActiveCamera()
		orientation = camera.GetOrientation()
		print('Initial orientation ', orientation)
		"""

	def delete_plots(self):
		for i in range( self.removable_actors.GetNumberOfItems() ):
			select_actor = self.removable_actors.GetItemAsObject(i)
			self.renderer.RemoveActor(select_actor)
		self.removable_actors.RemoveAllItems()

	def delete_billboard(self):
		for actor in self.style.text_Billboard_actors :
			self.renderer.RemoveActor(actor)
		self.style.text_Billboard_actors = []

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
			
			self.style.pickuptable_actors.AddItem(actor)
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

	def plot_text_2d(self, text = "2D text", pos = (10, 10), color = (1.0,1.0,1.0)):
		act_txt = vtk.vtkTextActor()
		act_txt.SetInput(text)
		prop_txt = act_txt.GetTextProperty()
		prop_txt.SetFontFamilyToArial()
		prop_txt.SetFontSize(18)
		prop_txt.SetColor(color)
		act_txt.SetDisplayPosition(*pos)
		self.renderer.AddActor(act_txt)

	def plot_text_Billboard(self, text = "3D text", pos = (1.0,1.0,1.0), color = (1.0,1.0,1.0), size = 0.005):
		act_txt = vtk.vtkBillboardTextActor3D()
		act_txt.SetPosition(pos)
		act_txt.SetInput(text)
		act_txt.SetScale(size)
		prop_txt = act_txt.GetTextProperty()
		prop_txt.SetFontFamilyToArial()
		prop_txt.SetFontSize(32)
		prop_txt.SetColor(color)
		self.renderer.AddActor(act_txt)
		self.style.text_Billboard_actors.append(act_txt)

	def plot_text_3d(self, text = "3D text", pos = (1.0,1.0,1.0), color = (1.0,1.0,1.0), size = 0.01):
		act_txt = vtk.vtkTextActor3D()
		act_txt.SetOrientation(90,90,0)# (0.0, -0.0, 0.0)
		act_txt.SetPosition(pos)
		act_txt.SetInput(text)
		act_txt.SetScale(size)
		prop_txt = act_txt.GetTextProperty()
		prop_txt.SetFontFamilyToArial()
		prop_txt.SetFontSize(32)
		prop_txt.SetColor(color)
		
		#prop_txt.SetBold(0)
		#prop_txt.SetItalic(0)
		#prop_txt.SetShadow(0)
		# 
		# prop_txt.SetVerticalJustificationToCentered()
		# 
		
		self.renderer.AddActor(act_txt)


	def plot_grid(self):
		axes = vtk.vtkAxesActor()
		axes_widget = vtk.vtkOrientationMarkerWidget()
		axes_widget.SetOutlineColor(1, 1, 1)
		axes_widget.SetOrientationMarker(axes)
		axes_widget.SetInteractor(inter)
		axes_widget.SetViewport(0., 0., 0.4, 0.4)
		axes_widget.EnabledOn()
		axes_widget.InteractiveOff()

	def plot_faces_dendrite(self, vertices, container_faces):

		c = plt.get_cmap('hsv', len(container_faces) ) 
		for i in range(len(container_faces)):
			actor = self.plot_mesh(vertices, container_faces[i], color=c(i)[:3] )
			self.removable_actors.AddItem(actor)
			self.renderer.AddActor(actor)


	def plot_spines(self):
		# Plot spine heads
		c = plt.get_cmap('summer', len( self.head_fs )+1 ) # hsv, flag
		for i, f in enumerate( self.head_fs ):
			if f is None:
				print('Null spine id: ', i)
				continue
			else:
				actor = self.plot_mesh(self.vertices, f, color = c(i)[:3] )
				self.renderer.AddActor(actor)


		# Plot spine necks
		c = plt.get_cmap('winter', len(self.neck_fs)+1 ) # hsv, flag
		for i, f in enumerate( self.neck_fs ):
			actor = self.plot_mesh(self.vertices, f, color = c(i)[:3] )
			self.renderer.AddActor(actor)


	def plot_disks(self, locations, normals, radiuses, color = (1.0, 1.0, 1.0) ):

		num = locations.shape[0]
		for i in range(num):
			self.plot_disk(locations[i,:], normals[i,:], radiuses[i], color = color)

	def plot_disk(self, location, normal, radius, color = (1.0, 1.0, 1.0) ):
		polygon_circle = vtk.vtkRegularPolygonSource()
		polygon_circle.SetNumberOfSides(10)
		polygon_circle.SetNormal( normal   )
		polygon_circle.SetRadius( radius   )
		polygon_circle.SetCenter( location )

		mapper_circle = vtk.vtkPolyDataMapper()
		mapper_circle.SetInputConnection(polygon_circle.GetOutputPort())

		actor = vtk.vtkActor()
		actor.SetMapper(mapper_circle)
		actor.GetProperty().SetColor(*color)
		#actor.GetProperty().SetLineWidth(width)
		self.renderer.AddActor(actor)



