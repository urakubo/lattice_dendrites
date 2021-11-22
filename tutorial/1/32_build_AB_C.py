import os
import numpy as np
import h5py
from pyLM.units import *
from pyLD import * 


input_morph_file = 'models/ball_and_stick.h5'
# Start variable
output_lm_file   = 'models/AB_C.lm'
# End variable

ext  = 'default'
cyt  = 'cytoplasm'
er   = 'er'
mito = 'mito'
domains = {ext: 0, cyt: 1, er: 2, mito: 3}


class SetMolecules:
	def __init__(self, cell):
		self.molecular_names  = ['A','B','C']
		self.cell = cell
		self.cell.define_species(self.molecular_names)

	def add_molecules(self, domain_name):
		concs = [1,100, 0]  # uM
		for name, conc in zip(self.molecular_names, concs):
			self.cell.add_molecule_uM(name, conc, domain_name )

	def set_diffusion(self, domain_name):
		diff_coeff = 10 * 1e-12 # Diffusion coefficient in SI
		for name in self.molecular_names:
			self.cell.set_diffusion(name, diff_coeff, domain_name)

	def set_reactions(self, domain_name):
		# AB binding to C
		kon = 1
		kof = 1
		elf.cell.reac_twoway_uM(reac=('A','B'), prod='C', rates=(kon, kof), domain_name=domain_name)


print('\nLoad geometry data.\n')
with h5py.File(input_morph_file,'r') as r:
	dendrite_not_mitochondrion_not_ER = r['dendrite not mitochondrion not ER'][()]
	pitch     = r['unit length (um)'][()]
	vol_PSD   = r['psd faces in volume'][()]
	vol_bound = r['bound faces in volume'][()]

volume = (dendrite_not_mitochondrion_not_ER > 0) * domains[cyt]
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