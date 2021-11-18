import os
import numpy as np
import h5py
from pyLM.units import *
from pyLD import * 


input_morph_file = 'models/ball_and_stick.h5'
output_lm_file   = 'models/photobleach.lm'

ext  = 'default'
cyt  = 'cytoplasm'
er   = 'er'
mito = 'mito'
domains = {ext: 0, cyt: 1, er: 2, mito: 3}


class SetMolecules:
	def __init__(self, cell):
		self.name_yfp = 'YFP'
		self.cell     = cell
		self.cell.define_species(self.name_yfp)

	def add_molecules(self, domain_name):
		conc_yfp = 1 # uM
		self.cell.add_molecule_uM(self.name_yfp, conc_yfp, domain_name)

	def set_diffusion(self, domain_name):
		d_yfp    = 1 * 1e-12 # (um2/s) Kang 2012; Traffic 13:1589-1600
		self.cell.set_diffusion(self.name_yfp, d_yfp, domain_name)

	def set_reactions(self, domain_name):
		pass


print('\nLoad geometry data.\n')
with h5py.File(input_morph_file,'r') as r:
	dendrite_not_mitochondrion_not_ER = r['dendrite not mitochondrion not ER'][()]
	pitch          = r['unit length per voxel (um)'][()]
	PSD            = r['PSD'][()]
	vol_bound      = r['boundary areas in volume'][()]

volume = (dendrite_not_mitochondrion_not_ER > 0) * domains[cyt]
vol_PSD  = vol_bound * (PSD > 0)
surfaces = {'PSD': vol_PSD, 'cell boundary': vol_bound}


print('\nBuild a lm model.\n')
spacing_in_m = micron(pitch) # Units in SI
cell = BuildAnyShape(volume, domains, spacing_in_m, surfaces)


print('\nSet molecules.\n')
molecules = SetMolecules(cell)
molecules.add_molecules(cyt)
molecules.set_diffusion(cyt)
molecules.set_reactions(cyt)


print('\nSet time.\n')
cell.sim.setTimestep(microsecond(3.0))
cell.sim.setWriteInterval(0.05)
cell.sim.setLatticeWriteInterval(0.05)
cell.sim.setSimulationTime(5.0)


print('\nSave simulation setup.\n')
if(os.path.isfile(output_lm_file) == True):
    os.remove(output_lm_file)
cell.sim.save(output_lm_file)