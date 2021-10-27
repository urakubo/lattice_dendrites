

from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
from skimage import morphology

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


def get_domain_concs(filenames, targs):
	"""Obtain time development of molecular concentrations.
	
	A 3D image is padded to a multiple of 32 voxels,
	because LM only accepts the size of volume.
	
	Args:
		filenames (list[str]): Output filenames of LM
		targs (list[str]): Target domains

	Returns:
		(tuple): Tuple containing:

		- time (numpy[float]): Time in s
		- concs (numpy[float]): Concentration in uM
		- numbers (numpy[int]): Numbers of Molecules
	"""

	for i, fname in enumerate(filenames):
	    with h5py.File(fname,'r') as f:
	        t  = np.array( f['t'][()] )
	        uM = f['conc in uM']
	        number = f['number']
	        molecules = list(uM.keys())
	        # Conc
	        tmp1 = np.zeros_like( uM[ molecules[0] ][()] )
	        tmp2 = np.zeros_like( number[ molecules[0] ][()] )
	        # print('tmp.shape: ', tmp.shape)
	        for targ in targs:
	            tmp1 += uM[targ][()]
	            tmp2 += number[targ][()]

	    # Connect
	    if i == 0:
	        #print('molecules: ', molecules)
	        concs   = tmp1
	        numbers = tmp2
	        #t[-1] = 10.0
	        time = t
	        #print('t: ', t)
	    else:
	        #print('uMs.shape : ', uMs.shape)
	        #print('tmp.shape: ', tmp.shape)
	        concs   = np.vstack( (concs, tmp1[1:,:]) )
	        numbers = np.vstack( (numbers, tmp2[1:,:]) )
	        time    = np.hstack( (time, t[1:]+Ts[-1]) )     
	    # print('No', i, ', Filename: ', fname)
	return time, concs, numbers


def get_species_name(filename):
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

		- num_voxels (int): Number of voxels of the target domain
		- volume_in_L (float): Volume of the target domain
		- spacing: Unit length (um)
		- s: Pairs of molecular name and id
	"""
	with h5py.File(filename,'r') as f:
	    data = f['Model']['Diffusion']['LatticeSites'][()]
	    spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
	    mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')

	## Volume
	if isinstance(domain_ids, int) | isinstance(domain_ids, bool) :
	    num_voxels  = np.count_nonzero(data == id_domains)
	    volume_in_L = num_voxels * spacing * spacing * spacing * 1000
	elif isinstance(domain_ids, list) | isinstance(domain_ids, tuple) :
	    num_voxels  = []
	    volume_in_L = []
	    for domain_id in domain_ids:
	        tmp_num = np.count_nonzero(data == domain_id)
	        tmp_vol = num_voxels * spacing * spacing * spacing * 1000
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


