

from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
from skimage import morphology
import h5py

def get_face_ids_inside(label_volume, v, f):
	"""Obtain the face ids that are located in label_volume.

	Args:
		label_volume (numpy[bool/int/float]): Target volume
		v[:,3] (float): Vertices
		f[:,3] (int)  : faces

	Returns:
		(id_face[int]): Ids of faces
	"""
	xvnum, yvnum, zvnum = label_volume.shape
	face_loc   = ( v[f[:,0]]+v[f[:,1]]+v[f[:,2]] ) / 3.0
	face_voxel = np.round( face_loc ).astype('int')
	face_voxel = (face_voxel < 0) + (face_voxel >= 0) * face_voxel
	face_voxel[:,0] = (face_voxel[:,0] >= xvnum) * (xvnum-1) + (face_voxel[:,0] < xvnum) * face_voxel[:,0]
	face_voxel[:,1] = (face_voxel[:,1] >= yvnum) * (yvnum-1) + (face_voxel[:,1] < yvnum) * face_voxel[:,1]
	face_voxel[:,2] = (face_voxel[:,2] >= zvnum) * (zvnum-1) + (face_voxel[:,2] < zvnum) * face_voxel[:,2]
	id_face = (label_volume[face_voxel[:,0],face_voxel[:,1],face_voxel[:,2]] > 0)

	return id_face

def get_volume(filename, id_domains):
	with h5py.File(filename,'r') as f:
		data = f['Model']['Diffusion']['LatticeSites'][()]
		Spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
	if isinstance(id_domains, int) | isinstance(id_domains, bool) :
		return (data == id_domains)
	else:
		print("id_domains must be int or bool.")
		return False

def smooth_volume_opening_then_closing(volume, radius_in_voxel=1):
	"""Smooth a volume by a one-round opening then closing.

	Args:
		volume (numpy[bool]): Target volume
		radius_in_voxel (int): Ball radius

	Returns:
		(numpy[bool]): Smoothing volume
	"""
	ball = morphology.ball(radius_in_voxel)
	volume = morphology.binary_opening(volume, ball) # morphology.ball(1)
	volume = morphology.binary_closing(volume, ball) # morphology.ball(1)
	return volume

def smooth_volume_closing_then_opening(volume, radius_in_voxel=1):
	"""Smooth a volume by a one-round closing then opening.

	Args:
		volume (numpy[bool]): Target volume
		radius_in_voxel (int): Ball radius

	Returns:
		(numpy[bool]): Smoothing volume
	"""
	ball = morphology.ball(radius_in_voxel)
	volume = morphology.binary_closing(volume, ball) # morphology.ball(1)
	volume = morphology.binary_opening(volume, ball) # morphology.ball(1)
	return volume

def smooth_volume_erosion_then_dilation(volume, radius_in_voxel=1):
	"""Smooth a volume by a one-round erosion then dilation.

	Args:
		volume (numpy[bool]): Target volume
		radius_in_voxel (int): Ball radius

	Returns:
		(numpy[bool]): Smoothing volume
	"""
	ball = morphology.ball(radius_in_voxel)
	volume = morphology.binary_erosion(volume, ball) # morphology.ball(1)
	volume = morphology.binary_dilation(volume, ball) # morphology.ball(1)
	return volume

def smooth_volume_dilation_then_erosion(volume, radius_in_voxel=1):
	"""Smooth a volume by a one-round dilation then erosion.

	Args:
		volume (numpy[bool]): Target volume
		radius_in_voxel (int): Ball radius

	Returns:
		(numpy[bool]): Smoothing volume
	"""
	ball = morphology.ball(radius_in_voxel)
	volume = morphology.binary_dilation(volume, ball) # morphology.ball(1)
	volume = morphology.binary_erosion(volume, ball) # morphology.ball(1)
	return volume

def lmpad(volume):
	"""3D padding into a multiple of 32 voxels, because LM only accepts the size of volume.

	Args:
		volume (numpy): Three dimentional numpy array

	Returns:
		(numpy): Padded volume
	"""

	nx,ny,nz = volume.shape
	min_lattice_size = 32 
	lx = np.ceil(1.0*nx/min_lattice_size)*min_lattice_size -nx
	ly = np.ceil(1.0*ny/min_lattice_size)*min_lattice_size -ny
	lz = np.ceil(1.0*nz/min_lattice_size)*min_lattice_size -nz
	lx1 = np.floor(lx/2)
	ly1 = np.floor(ly/2)
	lz1 = np.floor(lz/2)
	lx2 = lx - lx1
	ly2 = ly - ly1
	lz2 = lz - lz1
	padding = np.array( [[lx1,lx2],[ly1,ly2],[lz1,lz2]] , dtype='int')
	volume  = np.pad(volume, padding)
	return volume

def num_to_uM(num_molecules, num_voxels, spacing_in_m):
	"""Get concentration in uM from the number(s) of molecules.

	Args:
		num_molecules (int/float/numpy[int/float]): Number(s) of molecules
		num_voxels (int/float/numpy[int/float]): Number(s) of voxels.
		spacing (int/float): Spacing per lattice.
	Returns:
		(int/float/numpy[int/float]): Concentration in uM
	"""
	NA = 6.022e23
	number_per_1umol = NA /(1e6)
	volume_in_L = num_voxels * spacing_in_m * spacing_in_m * spacing_in_m * 1000
	number_per_1uM = number_per_1umol * volume_in_L
	conc = num_molecules / number_per_1uM
	return conc


def uM_to_num(conc_in_uM, num_voxels, spacing_in_m):
	"""Get the number(s) of molecules from the concentration(s) in uM.

	Args:
		conc_in_uM (int/float/numpy[int/float]): Conc(s) of molecules
		num_voxels (int/float/numpy[int/float]): Number(s) of voxels.
		spacing_in_m (int/float): Spacing per lattice.
	Returns:
		(int/numpy[int]): number of molecules
	"""
	NA = 6.022e23
	number_per_1umol = NA /(1e6)
	volume_in_L = num_voxels * spacing_in_m * spacing_in_m * spacing_in_m * 1000
	number_per_1uM = number_per_1umol * volume_in_L
	number = conc_in_uM * number_per_1uM
	number = np.int(np.floor(number))
	return number

def connect_total_concs(lm_filenames, species, domain_ids):
	"""Connect time developments of the number and concentration of specified molecules.
	Note that concentrations are simply obtained by the total numbers of molecules divided by the volume of a specified domain(s).

	Args:
		lm_filenames (str / list[str] / tuple[str]): Filename(s) of lm simulation.
		species (str / list[str] / tuple[str]): Molecular species. They are summed if multiple species are specified.
		domain_ids (int / list[int] / tuple[int]): Target domain ids. They are summed if multiple domains are specified.

	Returns:
		(tuple): Tuple containing:

		- timepoints (numpy[float]): Time in s
		- concs (numpy[float]): Concentration in uM
		- numbers (numpy[int]): Numbers of Molecules
	"""
	# Check arguments
	if isinstance(lm_filenames, str):
	    label_conc_filenames = [label_conc_filenames]
	elif isinstance(lm_filenames, list) | isinstance(lm_filenames, tuple) :
		pass
	else:
		raise ValueError('lm_filenames must be str, list, or tuple.')


	# Timepoints
	for i, fname in enumerate(lm_filenames):
	    t, c, n = get_total_concs(fname, species=species, domain_ids=domain_ids)
	    if i == 0:
	        timepoints = t
	        concs      = c
	        numbers    = n
	    else:
	        timepoints = np.hstack( (timepoints, t[1:]+timepoints[-1]) )
	        concs      = np.concatenate([concs   , c[1:] ], axis=0)
	        numbers    = np.concatenate([numbers , n[1:] ], axis=0)

	return timepoints, concs, numbers


def get_total_concs(lm_filename, species, domain_ids):
	"""Get time developments of the number and concentration of specified molecules.
	Note that concentrations are simply obtained by the total numbers of molecules divided by the volume of a specified domain(s).

	Args:
		lm_filename (str): Filename of lm simulation.
		species (str / list[str] / tuple[str]): Molecular species. They are summed if multiple species are specified.
		domain_ids (int / list[int] / tuple[int]): Target domain ids. They are summed if multiple domains are specified.

	Returns:
		(tuple): Tuple containing:

		- timepoints (numpy[float]): Time in s
		- concs (numpy[float]): Concentration in uM
		- numbers (numpy[int]): Numbers of Molecules
	"""

	if not isinstance(lm_filename, str):
		raise ValueError('lm_filename must be str.')

	if isinstance(species, str):
	    species = [species]
	elif isinstance(species, list) | isinstance(species, tuple) :
		pass
	else:
		raise ValueError('species must be str, list[str], or tuple[str].')

	if isinstance(domain_ids, int):
	    domain_ids = [domain_ids]
	elif isinstance(domain_ids, list) | isinstance(domain_ids, tuple) :
		pass
	else:
		raise ValueError('domain_ids must be int, list[int], or tuple[int].')

	with h5py.File(lm_filename, 'r') as file:
		timepoints = file['Simulations']['0000001']['LatticeTimes'][()]
		numbers = file['Simulations']['0000001']['SpeciesCounts'][()]
		spacing = file['Model']['Diffusion'].attrs['latticeSpacing']
		mnames  = file['Parameters'].attrs['speciesNames'].decode().split(',')
		volume  = file['Model']['Diffusion']['LatticeSites'][()]

	ids = [i for i, key in enumerate(mnames) if key in species ]
	numbers = np.sum(numbers[:, ids], axis = 1)

	num_voxels  = 0
	for domain_id in domain_ids:
		num_voxels += np.count_nonzero(volume == domain_id)

	concs = num_to_uM(numbers, num_voxels, spacing)

	return timepoints, concs, numbers


def get_spacing(filename):
	"""Retrieve spacing from a LM/LM-output file.

	Args:
		filename (str): Output filename of LM

	Returns:
		spacing (float): Unit length (um)
	"""
	with h5py.File(filename,'r') as f:
	    spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
	return spacing


def get_species_names(filename):
	"""Retrieve species names from a LM/LM-output file.

	Args:
		filename (str): Output filename of LM

	Returns:
		(cell[str] = id): Pairs of molecular name and id
	"""

	with h5py.File(filename,'r') as f:
	    mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
	s = {}
	for i in range(len(mnames)):
	    s[mnames[i]] = i+1
	return s


def get_volume_info(filename, domain_ids):
	"""Retrieve domain info from a LM/LM-output file.

	Args:
		filename (str): Output filename of LM
		domain_ids(int/list[int]/tuple[int]): Target domain ids

	Returns:
		(tuple): Tuple containing:

		- num_voxels (int/list[int]): Number of voxels of the target domain(s)
		- volume_in_L (float/list[float]): Volume of the target domain(s)
		- spacing (float): Unit length (um)
		- s: Pairs of molecular name and id
	"""
	with h5py.File(filename,'r') as f:
	    data = f['Model']['Diffusion']['LatticeSites'][()]
	    spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
	    mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')

	## Volume
	if isinstance(domain_ids, int) | isinstance(domain_ids, bool) :
	    num_voxels  = np.count_nonzero(data == domain_ids)
	    volume_in_L = num_voxels * spacing * spacing * spacing * 1000
	elif isinstance(domain_ids, list) | isinstance(domain_ids, tuple) :
	    num_voxels  = []
	    volume_in_L = []
	    for domain_id in domain_ids:
	        tmp_num = np.count_nonzero(data == domain_id)
	        tmp_vol = tmp_num * spacing * spacing * spacing * 1000
	        num_voxels.append(tmp_num)
	        volume_in_L.append(tmp_vol)
	    
	    ## Molecular names
	s = {}
	for i in range(len(mnames)):
	    s[mnames[i]] = i+1
	return num_voxels, volume_in_L, spacing, s


def get_annot_colors(filename, ids):
	"""Obtain colors of painted areas in the file of UNI-EM morphometric plugin.

	Args:
		filename (str): Filename of the paint npz file (the morphometric plugin).
		spine_ids(list[int]/tuple[int]): Target spine ids

	Returns:
		(list): List containing:

		- (R, G, B): colors (0-1 float) of target ids
	"""
	with open(filename,'rb') as f:
	    list = pickle.load(f)
	cols = []
	for id in ids:
	    c = [x for x in list['list'] if x['id'] == id]
	    r = c[0]['r']/256.0
	    g = c[0]['g']/256.0
	    b = c[0]['b']/256.0
	    cols.append((r,g,b))
	return cols



class Params():
	def __init__(self, user_path):

	    self.annotator_files_path      = user_path
	    self.volume_path               = os.path.join(self.annotator_files_path, 'volume')
	    self.volume_file               = os.path.join(self.volume_path, 'volume.hdf5')

	    self.skeletons_path            = os.path.join(self.annotator_files_path, 'skeletons')
	    self.surfaces_path             = os.path.join(self.annotator_files_path, 'surfaces')
	    self.skeletons_whole_path      = os.path.join(self.annotator_files_path, 'skeletons', 'whole')
	    self.surfaces_whole_path       = os.path.join(self.annotator_files_path, 'surfaces' , 'whole')
	    self.paint_path                = os.path.join(self.annotator_files_path, 'paint')

	    self.surfaces_segment_info_json_file 	    = os.path.join(self.surfaces_path, 'segmentInfo.json')
	    self.surfaces_volume_description_json_file 	= os.path.join(self.surfaces_path, 'VolumeDescription.json')


