import os
import numpy as np
import h5py
from tut1_functions import *
from pyLD import * 

input_model_file  = "models/photobleach.lm"
targ_name         = 'YFP'
color             = (1,1,0)

input_model_file  = "models/Ca_dynamics.lm"
targ_name         = 'inactive NMDAR'
#targ_name         = 'Pump'
color             = (1,0,0)

molecule  = get_species_names(input_model_file)
targ_id   = molecule[targ_name]

print('\nLoad molecules')
with h5py.File(input_model_file, 'r') as file:
	lattice = file['Model']['Diffusion']['Lattice'][:,:,:,:]
m = get_loc_molecule(lattice, targ_id)

print('\nShow molecules')
input_morpho_file = "models/ball_and_stick.h5"
mlab = show_dendrite(input_morpho_file)
mlab.points3d(m[0], m[1], m[2],\
	color = color,\
	scale_factor = 1.5)

mlab.savefig( os.path.join('imgs',targ_name+'.png') )
mlab.show()
