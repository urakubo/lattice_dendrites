

import sys, os, glob, pickle
import json
import h5py
import numpy as np
from skimage import morphology
import trimesh
import gzip
import pymeshfix
import pyvista as pv
from .utils import Params

main_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
##
class GenerateClosedVolumesFromUniEM():

	"""Generate volumes from the painted areas in UNI-EM annotator.

	Args:
	    folder_annot (str): UNI-EM annotator folder
	Returns:
		(pyLD. GenerateClosedVolumesFromUniEM):  GenerateClosedVolumesFromUniEM object
	"""
	def __init__(self, folder_annot):
		self.annot = Params(folder_annot)
		if self._check_annot_folder() == False:
			print('folder_annot is not a annotation folder: ', folder_annot)
			return False
		self._load_info()

	def _check_annot_folder(self):
		if  os.path.exists(self.annot.skeletons_whole_path) and \
		    os.path.exists(self.annot.surfaces_whole_path) and \
		    os.path.exists(self.annot.paint_path) and \
		    os.path.isfile(self.annot.surfaces_segment_info_json_file) and \
		    os.path.isfile(self.annot.surfaces_volume_description_json_file) and \
		    os.path.isfile(self.annot.volume_file) :
		    return True
		else:
		    return False

	def _load_info(self):
		#
		with open(self.annot.surfaces_volume_description_json_file, 'rb') as f:
			volume_info = json.load(f)
		# print(volume_info)
		self.ph = volume_info['pitch_um']['x']
		self.pw = volume_info['pitch_um']['y']
		self.pz = volume_info['pitch_um']['z']
		#
		self.whole_mesh_filenames    = glob.glob(os.path.join(self.annot.surfaces_whole_path, "*.stl"))
		self.whole_mesh_names_wo_ext = [ os.path.splitext(os.path.basename(text))[0] for text in self.whole_mesh_filenames ]
		print('Existent surface(s): ', ', '.join( self.whole_mesh_names_wo_ext ) )
		print()

	def generate(self, mesh_id, dilation_radius = 1):
		"""Generate volumes.

		Args:
		    mesh_id (int): Target mesh id
		    dilation_radius (int): dilation radius

		Returns:
			(tuple): Tuple containing:

			- new_labels (numpy[int]): Generated labels in voxel space (3D array) 
			- ids_new_labels (numpy[int]): Ids of new labels (1D array)
			- ref_volume (numpy[bool]): Referred id in voxel space (3D array)
		"""

		# Load volume
		with h5py.File( self.annot.volume_file, 'r') as f:	
			volume = f['volume'][()]
		ref_volume = ( volume == int(mesh_id) ).astype('bool')
		new_labels = np.zeros_like( ref_volume ).astype(np.int)

		# Generate filenames of surface meshes
		whole_mesh_filename     = os.path.join(self.annot.surfaces_whole_path, str(mesh_id).zfill(10)+".stl")
		whole_mesh_name_wo_ext  = os.path.splitext(os.path.basename(whole_mesh_filename))[0]

		# Check painted meshes
		part_mesh_name_wildcard = os.path.normpath(os.path.join(self.annot.paint_path, whole_mesh_name_wo_ext+"-*.pickle"))

		# print('whole_mesh_filename    : ', whole_mesh_filename)
		# print('whole_mesh_name_wo_ext : ', whole_mesh_name_wo_ext)
		# print('part_mesh_name_wildcard: ', part_mesh_name_wildcard)

		part_mesh_filenames = glob.glob(part_mesh_name_wildcard)
		if part_mesh_filenames == [] :
			print('No paint areas found.')
			return False
		part_meshes_wo_ext = [ os.path.splitext(os.path.basename(text))[0] for text in part_mesh_filenames ]
		part_mesh_ids      = [ int( text.split('-')[1] ) for text in part_meshes_wo_ext ]

		# Load surface meshes
		whole_mesh = trimesh.load( whole_mesh_filename )
		whole_mesh.vertices[:,0] /= self.pw
		whole_mesh.vertices[:,1] /= self.ph
		whole_mesh.vertices[:,2] /= self.pz

		ids_new_labels = []
		for part_mesh_filename, part_mesh_id in zip(part_mesh_filenames, part_mesh_ids) :
			with open(part_mesh_filename, 'rb') as file:
				data = pickle.load(file)
			closed_mesh = self._get_closed_mesh(whole_mesh, data['painted'])
			if closed_mesh.volume is None :
				continue
			#
			ids_new_labels.append(part_mesh_id)
			v = self._get_volume(closed_mesh)

			# Embed it into a whole gridspace
			wmin = 1 + np.floor( np.min( closed_mesh.vertices[:,0] )).astype(int)
			hmin = 1 + np.floor( np.min( closed_mesh.vertices[:,1] )).astype(int)
			zmin = 1 + np.floor( np.min( closed_mesh.vertices[:,2] )).astype(int)

			wnum = v.shape[0]
			hnum = v.shape[1]
			znum = v.shape[2]

			print('ID: ', part_mesh_id,' Volume:', closed_mesh.volume)
			print('wnum, hnum, znum : ', wnum, hnum, znum)
			print('Volume v: ', np.sum(v.astype('int')))

			tmp_labels = np.zeros_like( ref_volume, dtype='bool' )
			tmp_labels[wmin:wmin+wnum , hmin:hmin+hnum, zmin:zmin+znum] = v.astype('bool')
			tmp_labels = morphology.binary_dilation(tmp_labels, morphology.ball( dilation_radius, dtype='bool' ) )
			new_labels[tmp_labels & ref_volume] = part_mesh_id
			# new_labels[tmp_labels > 0] = part_mesh_id

		return new_labels, ids_new_labels, ref_volume



	def _get_closed_mesh(self, mesh, data):

		f = mesh.faces
		unzipped_tri = gzip.decompress(data)
		sub_face_id = []
		for i in range( f.shape[0] ) :
			if (unzipped_tri[i*3:i*3+3] == b'\x01\x01\x01') :
				sub_face_id.append(f[i,:])

		mesh = trimesh.Trimesh(mesh.vertices, np.array(sub_face_id))
		mesh.merge_vertices()
		mesh.remove_degenerate_faces()
		mesh.remove_duplicate_faces()
		mesh.remove_unreferenced_vertices()

		vclean = mesh.vertices
		fclean = mesh.faces
		part_mesh = pymeshfix.MeshFix(vclean, fclean)
		part_mesh.repair()
		part_mesh.plot() # Visualization of cloased meshes

		closed_mesh = trimesh.Trimesh(vertices=part_mesh.v, faces=part_mesh.f)
		closed_mesh.merge_vertices()
		closed_mesh.remove_degenerate_faces()
		closed_mesh.remove_duplicate_faces()
		closed_mesh.remove_unreferenced_vertices()
		return closed_mesh


	def _get_volume(self, mesh):

		verts = mesh.vertices
		faces = mesh.faces

		wmin = np.floor( np.min(verts[:,0])).astype(int)
		hmin = np.floor( np.min(verts[:,1])).astype(int)
		zmin = np.floor( np.min(verts[:,2])).astype(int)

		wmax = np.floor( np.max(verts[:,0])).astype(int)
		hmax = np.floor( np.max(verts[:,1])).astype(int)
		zmax = np.floor( np.max(verts[:,2])).astype(int)

		num = faces.shape[0]
		faces = np.hstack([np.ones([num,1]).astype(int)*3,faces])
		surf = pv.PolyData(verts, faces)
		ix = np.arange(wmin, wmax, 1)
		iy = np.arange(hmin, hmax, 1)
		iz = np.arange(zmin, zmax, 1)

		x, y, z = np.meshgrid(ix, iy, iz)
		grid = pv.StructuredGrid(x, y, z)
		ugrid = pv.UnstructuredGrid(grid)
		selection = ugrid.select_enclosed_points(surf, tolerance=0.0, check_surface=False)
		mask = selection.point_arrays['SelectedPoints'].view(np.bool)
		voxels = mask.reshape([iz.shape[0] , ix.shape[0], iy.shape[0] ])
		voxels = voxels.transpose((1, 2, 0))
		return voxels




