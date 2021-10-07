
import numpy as np


def lmpad(volume):
	"""
	A 3D image is padded to a multiple of 32 voxels,
	because LM only accepts the size of volume.
	
	Args:
		volume (any):	3 dimentional numpy array
	
	Returns:
		volume : padded volume
	
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


def get_domain_concs(filenames, Targs):
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
            for Targ in Targs:
                tmp1 += uM[Targ][()]
                tmp2 += number[Targ][()]

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
    with h5py.File(filename,'r') as f:
        mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
    S = {}
    for i in range(len(mnames)):
        S[mnames[i]] = i+1
    return S


def get_volume_info(filename, id_domains):
    with h5py.File(filename,'r') as f:
        data = f['Model']['Diffusion']['LatticeSites'][()]
        Spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
        mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')

    ## Volume
    if isinstance(id_domains, int) | isinstance(id_domains, bool) :
        num_voxels  = np.count_nonzero(data == id_domains)
        volume_in_L = num_voxels * Spacing * Spacing * Spacing * 1000
    elif isinstance(id_domains, list) | isinstance(id_domains, tuple) :
        num_voxels  = []
        volume_in_L = []
        for id_domain in id_domains:
            tmp_num = np.count_nonzero(data == id_domain)
            tmp_vol = num_voxels * Spacing * Spacing * Spacing * 1000
            num_voxels.append(tmp_num)
            volume_in_L.append(tmp_vol)
        
        ## Molecular names
    S = {}
    for i in range(len(mnames)):
        S[mnames[i]] = i+1
    return num_voxels, volume_in_L, Spacing, S

    
def get_annot_colors(filename, ids_spine):
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


