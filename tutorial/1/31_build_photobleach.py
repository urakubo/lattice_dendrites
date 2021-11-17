import os
import numpy as np
import h5py
from pyLM import *
from pyLM.units import *
from pyLD import * 


input_filename_morph = 'models/ball_and_stick.h5'
output_filename_lm   = 'models/ball_and_stick_photobleach.lm'

ext  = 'default'
cyt  = 'cytoplasm'
er   = 'er'
mito = 'mito'
domains = {ext: 0, cyt: 1, er: 2, mito: 3}


print('\nLoad geometry data.\n')

with h5py.File(input_filename_morph,'r') as r:
    dendrite_not_mitochondrion_not_ER = r['dendrite not mitochondrion not ER'][()]
    pitch          = r['unit length per voxel (um)'][()]
    PSD            = r['PSD'][()]
    vol_bound      = r['boundary areas in volume'][()]

volume = (dendrite_not_mitochondrion_not_ER > 0) * domains[cyt]
vol_PSD  = vol_bound * (PSD > 0)
surfaces = {'PSD': vol_PSD, 'cell boundary': vol_bound}


print('\nBuild a lm model.\n')
spacing=micron(pitch) # Units in SI
nx,ny,nz = volume.shape
sim  = RDME.RDMESimulation(dimensions=micron(nx*pitch,ny*pitch,nz*pitch), spacing=spacing)
cell = BuildAnyShape(sim, volume, domains, surfaces)

num_voxel_cyt  = np.count_nonzero(volume == domains[cyt])
number_per_1uM = uM_to_num(conc_in_uM = 1, num_voxels = num_voxel_cyt, spacing_in_m = spacing)
print('Number of molecules per 1 uM:', number_per_1uM)


print('\nSet molecules.\n')
name_yfp = 'YFP'
d_YFP    = 1 * 1e-12 ## (um2/s) Kang 2012; Traffic 13:1589-1600
conc_YFP = 1 # uM
num_YFP  = conc_YFP * number_per_1uM
print('num_YFP : ', num_YFP)

sim.defineSpecies( [name_yfp] )
sim.modifyRegion( cyt ).setDiffusionRate( name_yfp , rate=d_YFP)

cell.add_solute_molecules( name_yfp, num_YFP, cyt) # Absolute number

print('\nSet time.\n')
sim.setTimestep(microsecond(3.0))
sim.setWriteInterval(0.05)
sim.setLatticeWriteInterval(0.05)
sim.setSimulationTime(5.0)

print('\nSave simulation setup.\n')
if(os.path.isfile(output_filename_lm) == True):
    os.remove(output_filename_lm)
sim.save(output_filename_lm)