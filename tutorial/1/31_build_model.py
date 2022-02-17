import os
import numpy as np
import h5py
from pyLM.units import micron, microsecond
from pyLD import * 

# FRAP
from set_molecule_FRAP import set_molecules
output_lm_file   = 'models/photobleach.lm'

'''
# Ca2+ influx via NMDARs
from set_molecule_Ca import set_molecules
output_lm_file   = 'models/Ca_dynamics.lm'
'''

print('\nLoad geometry data.')
input_morph_file = 'models/ball_and_stick.h5'
with h5py.File(input_morph_file,'r') as r:
	pitch_in_um = r['unit length (um)'][()]
	vol_cytosol = r['dendrite not mitochondrion not ER'][()]
	vol_PSD     = r['psd faces in volume'][()]
	vol_bound   = r['bound faces in volume'][()]

vol_cytosol = (vol_cytosol > 0).astype('uint8')
domains     = {'not cytosol': 0, 'cytosol': 1}
surfaces    = {'PSD': vol_PSD, 'cell boundary': vol_bound}


print('\nBuild a lm model.')
pitch_in_m = micron(pitch_in_um) # Unit in SI
cell = BuildAnyShape(vol_cytosol, domains, pitch_in_m, surfaces)
molecules = set_molecules(cell)


print('\nSet time.')
cell.sim.setTimestep(microsecond(3.0))
cell.sim.setWriteInterval(0.05)
cell.sim.setLatticeWriteInterval(0.05)
cell.sim.setSimulationTime(5.0)


print('\nSave overall setup.')
if(os.path.isfile(output_lm_file) == True):
    os.remove(output_lm_file)
cell.sim.save(output_lm_file)
