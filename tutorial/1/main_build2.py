import os
import numpy as np
import h5py
from pyLM import *
from pyLM.units import *
from pyLD import * 


input_filename_morph = 'models/ball_and_stick.h5'
output_filename_lm   = 'models/ball_and_stick_AB_C.lm'

ext  = 'default'
cyt  = 'cytoplasm'
er   = 'er'
mito = 'mito'
domains = {ext: 0, cyt: 1, er: 2, mito: 3}



class setMolecules:
    def __init__(self, cell):
        self.cell     = cell
        self.molcules  = ['A','B','C']
        self.cell.define_species(self.molcules)

    def addMolecules(self, domain_name):
        conc_A = 1   # uM
        conc_B = 100 # uM
        conc_C = 0
        cell.add_solute_molecules_uM( 'A', conc_A, domain_name ) # Absolute number
        cell.add_solute_molecules_uM( 'B', conc_B, domain_name ) # Absolute number
        cell.add_solute_molecules_uM( 'C', conc_C, domain_name ) # Absolute number

    def setDiffusion(self, domain_name):
        d_ABC = 10 * 1e-12 # Diffusion constant
        self.d_mol = {'A':d_ABC ,
                      'B':d_ABC ,
                      'C':d_ABC }
        for k, v in self.d_mol.items():
            self.cell.set_diffusion_rate(k, v, domain_name)

    def setReactions(self, domain_name):
        # uM_s_to_number_s
        on = self.cell.per_uM

        # Ca-CaM binding reactions
        kon_AB_C = 1 * on / 100
        kof_AB_C = 1 * 1

        # AB binding to C
        A = 'A'
        B = 'B'
        C = 'C'
        cyt = self.cell.modifyRegion( domain_name )
        cyt.addReaction(reactant=(A,B), product=C    , rate=kon_AB_C)
        cyt.addReaction(reactant= C   , product=(A,B), rate=kof_AB_C)



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


print('\nSet molecules.\n')
molecules = setMolecules(cell)
molecules.addMolecules(cyt)
molecules.setDiffusion(cyt)
molecules.setReactions(cyt)


print('\nSet time.\n')
sim.setTimestep(microsecond(3.0))
sim.setWriteInterval(0.05)
sim.setLatticeWriteInterval(0.05)
sim.setSimulationTime(5.0)

print('\nSave simulation setup.\n')
if(os.path.isfile(output_filename_lm) == True):
    os.remove(output_filename_lm)
sim.save(output_filename_lm)
