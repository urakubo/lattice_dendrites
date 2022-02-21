import sys, os
import h5py
from pyLD import save_uniem_annotator

filename     = 'models/realistic_dendrite.h5'
annot_folder = 'annot_realistic'

with h5py.File(filename,'r') as f:
	bound_v   = f['bound vertices'][()]
	bound_f   = f['bound faces'][()]
	mito_v    = f['mito vertices'][()]
	mito_f    = f['mito faces'][()]
	er_v      = f['er vertices'][()]
	er_f      = f['er faces'][()]
	volume    = f['volume'][()]
	pitch     = f['unit length (um)'][()]


bound_color = [192,192,192]
mito_color  = [255,255,152]
er_color    = [179,255,179]
surfaces = {1: [bound_v, bound_f, bound_color],\
			2: [mito_v , mito_f , mito_color],\
			3: [er_v   , er_f   , er_color]}

save_uniem_annotator(annot_folder, pitch, volume.astype('uint16'), surfaces)