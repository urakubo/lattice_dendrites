from __future__ import print_function
from __future__ import division

import h5py
import numpy as np
import os, sys, shutil
import subprocess as s
from os.path import join 
from lib.Utils import normalize_volume_size, get_species_name
gpu_id = '7'
targ_dir                = 'lms_3_4_2'
filename_lm             = join(targ_dir, '_model.lm')
filename_prerun         =[join(targ_dir, '_prerun_result.lm')]
filename_stimrun_prefix = join(targ_dir, '_stimrun_')
input_label_file = "lm_annot/labels.hdf5"

if not os.path.isfile(filename_lm):
    print('No model file         : ', filename_lm)
    sys.exit()
if not os.path.isfile(filename_prerun[0]):
    print('No prerun file     : ', filename_prerun)


## Load species name
S = get_species_name(filename_prerun[-1])

## Load spine labels    
with h5py.File(input_label_file,'r') as f:
    label_spine = f['dendrite'][()]

label_spine = normalize_volume_size( label_spine )
ids_spine, nums_spine_voxels = np.unique(label_spine, return_counts=True)
ids_spine         = ids_spine[1:-2]
#nums_spine_voxels = nums_spine_voxels[1:-2]
tmp = ids_spine[1::4]
ids_targ_spines = list( set(ids_spine) - set(tmp) )

print("spines id     : ", ids_spine)
print("targ_spines id: ", ids_targ_spines)
mask_spine = np.zeros_like(label_spine, dtype=np.bool)
for id_targ_spine in ids_targ_spines:
    mask_spine = (mask_spine | (label_spine == id_targ_spine))

# sys.exit()

#period = np.array('1.0  ',dtype=str)
period = np.string_(['1.0'])
#for i in range(10):
for i in range(30):

    with h5py.File(filename_prerun[-1],'r') as f:
        TimePoints = list(f['Simulations']['0000001']['Lattice'].keys())
        TimePoints.sort()
        print('Time ID: ', TimePoints[-1])
        Lattice4 = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
        SpeciesCount = f['Simulations']['0000001']['SpeciesCounts'][-1,:]
        num_targ_NR = 0
        for j in range(Lattice4.shape[3]):
            Lattice3 = Lattice4[:,:,:,j]
            num_targ_NR += np.count_nonzero( Lattice3[(Lattice3 == S['NR'] ) & mask_spine] )
            Lattice3[(Lattice3 == S['NR']) & mask_spine] = S['NR_Glu']
            Lattice4[:,:,:,j] = Lattice3 # Unnecessary?
            
        SpeciesCount[S['NR_Glu']-1] += num_targ_NR
        SpeciesCount[S['NR']-1]     -= num_targ_NR
        id = '%02d' % i
        filename = filename_stimrun_prefix + id +'.lm'
        filename_prerun.append(filename)
        print('Create stimrun file: ', filename)
        shutil.copy(filename_lm, filename)

        with h5py.File(filename,'a') as g:
            g['Model']['Diffusion']['Lattice'][()] = Lattice4
            g['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount
            g['Parameters'].attrs['maxTime'] = period
            g['Model']['Reaction']['ReactionRateConstants'][78,0] = 6000.0
        com = ['lm','-nr', '-g', gpu_id, '-r', '1', '-cr','5', '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename]
        print(' '.join(com))
        s.call(com)
