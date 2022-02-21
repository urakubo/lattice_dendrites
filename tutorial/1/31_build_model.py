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
	volume      = r['volume'][()]
	face_PSD    = r['psd faces in volume'][()]
	face_bound  = r['bound faces in volume'][()]

domains     = {'ecs': 0, 'cytosol': 1, 'mito': 2, 'er': 3}
surfaces    = {'PSD': face_PSD, 'cell boundary': face_bound}


print('\nBuild a lm model.')
pitch_in_m = micron(pitch_in_um) # Unit in SI
cell = BuildAnyShape(cytosol, domains, pitch_in_m, surfaces)
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
