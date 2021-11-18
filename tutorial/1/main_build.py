import os
import numpy as np
import h5py
from pyLM.units import *
from pyLD import * 


input_morph_file = 'models/ball_and_stick.h5'
output_lm_file   = 'models/AB_C.lm'

ext  = 'default'
cyt  = 'cytoplasm'
er   = 'er'
mito = 'mito'
domains = {ext: 0, cyt: 1, er: 2, mito: 3}


class SetMolecules:
    def __init__(self, cell):
        molecules  = ['A','B','C']
        self.cell  = cell
        self.cell.define_species(molecules)

    def add_molecules(self, domain_name):
        number_per_1uM = self.cell.number_per_1uM(domain_name)
        num_A = 1 * number_per_1uM
        num_B = 1 * number_per_1uM
        num_C = 0 * number_per_1uM
        self.cell.add_molecule( 'A', num_A, domain_name ) # Absolute number
        self.cell.add_molecule( 'B', num_B, domain_name ) # Absolute number
        self.cell.add_molecule( 'C', num_C, domain_name ) # Absolute number

    def set_diffusion(self, domain_name):
        d_ABC = 10  * 1e-12 # Diffusion constant
        d_mol = {'A':d_ABC ,
                 'B':d_ABC ,
                 'C':d_ABC }
        for k, v in d_mol.items():
            self.cell.set_diffusion_rate(k, v, domain_name)

    def set_reactions(self, domain_name):

        # uM_s_to_number_s
        on = self.cell.per_uM

        # Ca-CaM binding reactions
        kon_AB_C = 1 * on
        kof_AB_C = 1 * 1

        # AB binding to C
        A = 'A'
        B = 'B'
        C = 'C'
        cyt = self.cell.modifyRegion( domain_name )
        cyt.addReaction(reactant=(A,B), product=C    , rate=kon_AB_C)
        cyt.addReaction(reactant= C   , product=(A,B), rate=kof_AB_C)


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
cell = BuildAnyShape(sim, volume, domains, spacing_in_m, surfaces)


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