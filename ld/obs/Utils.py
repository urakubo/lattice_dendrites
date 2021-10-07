from __future__ import print_function
from __future__ import division

import numpy as np
import os
import pickle
import h5py
from os.path import join 

def get_active_synapse(name):
    ids_spine = list(range(1,23))
    if  name in ['lms_1_4']:
        ids_targ_spines = ids_spine[::4]
        ids_other_spines = list( set(ids_spine) - set(ids_targ_spines) )
    elif name in ['lms_2_4','lms300_2_4','lms450_2_4','lms600_2_4','lms900_2_4','lms000_2_4','lms1200_2_4']:
        ids_targ_spines = ids_spine[::2]
        ids_other_spines = list( set(ids_spine) - set(ids_targ_spines) )
    elif name in ['lms_3_4'] :
        ids_other_spines = ids_spine[::4]
        ids_targ_spines = list( set(ids_spine) - set(ids_other_spines) )
    elif name in ['lms_1_4_2'] :
        ids_targ_spines = ids_spine[1::4]
        ids_other_spines = list( set(ids_spine) - set(ids_targ_spines) )
    elif name in ['lms_2_4_2','lms300_2_4_2','lms450_2_4_2','lms600_2_4_2','lms900_2_4_2','lms000_2_4_2','lms1200_2_4_2'] :
        ids_targ_spines = ids_spine[1::2]
        ids_other_spines = list( set(ids_spine) - set(ids_targ_spines) )
    elif name in ['lms_3_4_2']:
        ids_other_spines = ids_spine[1::4]
        ids_targ_spines = list( set(ids_spine) - set(ids_other_spines) )
    else :
        print("get_active_synapse: no target name, ", name)
        return False
    return ids_spine, ids_targ_spines, ids_other_spines


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


def normalize_volume_size(_volume):
    nx,ny,nz = _volume.shape
    min_lattice_size = 32
    lx = np.ceil(1.0*nx/min_lattice_size)*min_lattice_size
    ly = np.ceil(1.0*ny/min_lattice_size)*min_lattice_size
    lz = np.ceil(1.0*nz/min_lattice_size)*min_lattice_size
    lx = lx.astype(np.int)
    ly = ly.astype(np.int)
    lz = lz.astype(np.int)
    volume = np.zeros((lx,ly,lz), dtype=_volume.dtype)
    volume[:nx,:ny,:nz] = _volume
    return volume

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
