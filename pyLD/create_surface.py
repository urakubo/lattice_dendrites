from __future__ import print_function
from __future__ import division

from skimage import measure
from skimage import morphology
import trimesh
import h5py
import numpy as np
import os

def create_surface(pitch, volume, PSD = None, num_smoothing = 30, method_smoothing = 'humphrey'):

	"""Create a smooth surface mesh of a volume.
	Smoothing is based a humphrey filter, and surface area per face is calcudated.
	Faces within PSD regions are labeled (Note, currently BEFORE smoothing).

	Args:
	    pitch (float): Unit length per voxel
	    volume (numpy): Input volume (3D array, bool preferred)
	    PSD (numpy): PSD volume  (3D array)
	    num_smoothing (int): Number of smoothing rounds for the target surface mesh
		method_smoothing(str): Smoothing method: 'laplacian', 'humphrey', 'mut_dif_laplacian', or 'taubin'

	Returns:
		(tuple): Tuple containing:

		- smooth_vertices (numpy[float]): Vertices of smoothing mesh (3xX array)
		- smooth_faces (numpy[int]): Faces of smoothing mesh (3xY array)
		- smooth_area_per_face (numpy[float]): Areas of faces (Y array)
		- surface_areas (numpy[float]): Surface areas in voxel space (3D array)
		- id_face_psd (numpy[bool]): Faces that are located at PSD (X array, bool)
	"""


	volume = volume.astype(np.bool)
	xvnum, yvnum, zvnum = volume.shape

	# Obtain the region of membrane
	cytosol = morphology.binary_erosion(volume, morphology.ball(1))
	memb    = volume ^ cytosol

	# Generate linked list
	memb_id = np.array(np.where(memb)).T

	# Marching cube
	v_march, f_march, normals, values = measure.marching_cubes(volume, 0.5, spacing=(1,1,1))
	v_march = v_march - 1


	# Smoothing
	trimesh.constants.tol.merge = 1e-7
	mesh = trimesh.Trimesh(vertices=v_march, faces=f_march)
	v1  = mesh.vertices
	v1center = np.sum(v1,0) / v1.shape[0]

	if method_smoothing == 'laplacian':
		mesh_smooth = trimesh.smoothing.filter_laplacian(mesh, iterations=num_smoothing)
		v2 = mesh_smooth.vertices
		v2center = np.sum(v2,0) / v2.shape[0]
		mesh_smooth.vertices += v1center - v2center
	elif method_smoothing == 'humphrey':
		mesh_smooth = trimesh.smoothing.filter_humphrey(mesh, alpha = 1.0, beta=0.0, iterations=num_smoothing)
	elif method_smoothing == 'mut_dif_laplacian':
		mesh_smooth = trimesh.smoothing.filter_mut_dif_laplacian(mesh, iterations=num_smoothing)
	elif method_smoothing == 'taubin':
		mesh_smooth = trimesh.smoothing.filter_taubin(mesh, iterations=num_smoothing)
	else :
		print('No smoothing method: ', method_smoothing)
		return False


	mesh_smooth.merge_vertices()
	mesh_smooth.remove_degenerate_faces()
	mesh_smooth.remove_duplicate_faces()
	v_smooth = mesh_smooth.vertices
	f_smooth = mesh_smooth.faces

	# PSD
	if PSD is not None:
		xvnum, yvnum, zvnum = PSD.shape
		face_loc   = ( v_smooth[f_smooth[:,0]]+v_smooth[f_smooth[:,1]]+v_smooth[f_smooth[:,2]] ) / 3.0
		face_voxel = np.round( face_loc ).astype(np.int)
		face_voxel = (face_voxel < 0) + (face_voxel >= 0) * face_voxel
		face_voxel[:,0] = (face_voxel[:,0] >= xvnum) * (xvnum-1) + (face_voxel[:,0] < xvnum) * face_voxel[:,0]
		face_voxel[:,1] = (face_voxel[:,1] >= yvnum) * (yvnum-1) + (face_voxel[:,1] < yvnum) * face_voxel[:,1]
		face_voxel[:,2] = (face_voxel[:,2] >= zvnum) * (zvnum-1) + (face_voxel[:,2] < zvnum) * face_voxel[:,2]
		id_face_psd = (PSD[face_voxel[:,0],face_voxel[:,1],face_voxel[:,2]] > 0)


	# Obrain the surface area of each trianglar face.
	f_vec_ac = v_smooth[f_smooth[:,0]] - v_smooth[f_smooth[:,2]]
	f_vec_bc = v_smooth[f_smooth[:,1]] - v_smooth[f_smooth[:,2]]
	f_outer  = np.cross(f_vec_ac, f_vec_bc) 
	f_areas  = np.linalg.norm(f_outer, axis=1) / 2
	# np.sum(f_areas) == measure.mesh_surface_area(v_smooth,f_smooth)

	# Obtain the location ids of trianglar faces in the voxel space
	face_loc  = ( v_smooth[f_smooth[:,0]]+v_smooth[f_smooth[:,1]]+v_smooth[f_smooth[:,2]] ) / 3.0
	face_voxel = np.round( face_loc ).astype(np.int)
	face_voxel = (face_voxel < 0) + (face_voxel >= 0) * face_voxel
	face_voxel[:,0] = (face_voxel[:,0] >= xvnum) * (xvnum-1) + (face_voxel[:,0] < xvnum) * face_voxel[:,0]
	face_voxel[:,1] = (face_voxel[:,1] >= yvnum) * (yvnum-1) + (face_voxel[:,1] < yvnum) * face_voxel[:,1]
	face_voxel[:,2] = (face_voxel[:,2] >= zvnum) * (zvnum-1) + (face_voxel[:,2] < zvnum) * face_voxel[:,2]

	# Set membrane area in the voxel space
	memb_area = np.zeros_like(volume, dtype=np.float)
	for i in range(f_areas.shape[0]):
		j = 1
		while (j < 20):
			loc = np.where(morphology.ball(j, dtype=np.bool))
			loc = np.array(loc).T - j + face_voxel[i,:]
			loc = loc[np.all(loc >= 0, axis=1),:]
			loc = loc[loc[:,0] < xvnum,:]
			loc = loc[loc[:,1] < yvnum,:]
			loc = loc[loc[:,2] < zvnum,:]
			tmp = memb[loc[:,0],loc[:,1],loc[:,2]]
			if np.sum(tmp) > 0:
				memb_area[loc[:,0],loc[:,1],loc[:,2]] += tmp.astype(np.float) / np.sum(tmp) * f_areas[i] # f_areas[i]
				break
			else :
				j = j + 1
				# print(j, i)
				# print(tmp)
				# print(loc)

	surface_areas  = memb_area * pitch * pitch
	smooth_vertices = v_smooth * pitch
	smooth_faces    = f_smooth
	smooth_area_per_face = f_areas * pitch * pitch
	if PSD is not None:
		return smooth_vertices, smooth_faces, smooth_area_per_face, surface_areas, id_face_psd
	else:
		return smooth_vertices, smooth_faces, smooth_area_per_face, surface_areas

