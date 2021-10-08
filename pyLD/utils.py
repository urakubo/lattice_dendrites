
import numpy as np


def lmpad(volume):
	"""3D padding for lattice microbes.

	A 3D image is padded to a multiple of 32 voxels,
	because LM only accepts the size of volume.

	Args:
		volume (any): Three dimentional numpy array

	Returns:
		volume: Padded volume

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
	padding = np.array( [[lx1,lx2],[ly1,ly2],[lz1,lz2]] , dtype=volume.dtype)
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
		Ts (numpy[float]): Time in s?, 
		uMs (numpy[float]): Concentration in uM.
		numbers (numpy[int]): Molecular number
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
	        uMs     = tmp1
	        numbers = tmp2
	        #t[-1] = 10.0
	        Ts  = t
	        #print('t: ', t)
	    else:
	        #print('uMs.shape : ', uMs.shape)
	        #print('tmp.shape: ', tmp.shape)
	        uMs     = np.vstack( (uMs, tmp1[1:,:]) )
	        numbers = np.vstack( (numbers, tmp2[1:,:]) )
	        Ts      = np.hstack( (Ts, t[1:]+Ts[-1]) )     
	    # print('No', i, ', Filename: ', fname)
	return Ts, uMs, numbers


def get_species_name(filename):
	"""Retrieve species names from a LM output file.

	Args:
		filename (str): Output filename of LM

	Returns:
		s (cell[str] = id): Molecular names and id
	"""

	with h5py.File(filename,'r') as f:
	    mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
	s = {}
	for i in range(len(mnames)):
	    s[mnames[i]] = i+1
	return s


def get_volume_info(filename, domain_ids):
	"""Retrieve domain info from a LM output file.

	Args:
		filename (str): Output filename of LM
		domain_ids(int/list[int]/tuple[int]): Target domain ids

	Returns:
		num_voxels (int): Number of voxels of the target domain.
		volume_in_L (float): Volume of the target domain.
		spacing: Unit length (unit?)
		s: Molecular names that specifies id
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


def get_annot_colors(filename, spine_ids):
	"""Obtain spine color from the UNI-EM morphometric plugin.

	Args:
		filename (str): Filename of the paint npz file (the morphometric plugin).
		spine_ids(list[int]/tuple[int]): Target spine ids

	Returns:
		cols (list[tuple[int]]): List of the RGB colors of the target ids.
	"""
	with open(filename,'rb') as f:
	    list = pickle.load(f)
	cols = []
	for id_spine in ids_spine:
	    c = [x for x in list['list'] if x['id'] == id_spine]
	    r = c[0]['r']/256.0
	    g = c[0]['g']/256.0
	    b = c[0]['b']/256.0
	    cols.append((r,g,b))
	return cols


