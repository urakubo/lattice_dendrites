import sys, os
import numpy as np
import h5py
from skimage import morphology
from pyLD import *
from tut1_functions import make_cylinder, add_shape


# Spiny dendrite
spine_head = morphology.ball(radius = 12)
spine_neck = make_cylinder(radius = 5, length = 50, direction = 2)
dendrite = make_cylinder(radius = 25, length = 90, direction = 0)

vol_dend = np.zeros((96,60,96), dtype=np.uint8)
vol_dend = add_shape(vol_dend, spine_head, [48,30,76])
vol_dend = add_shape(vol_dend, spine_neck, [48,30,51])
vol_dend = add_shape(vol_dend, dendrite  , [48,30,26])


# PSD
psd = morphology.ball(radius = 7)
vol_psd = add_shape(np.zeros_like(vol_dend), psd, [48,30,88])


# Mito
mito = make_cylinder(radius = 10, length = 90, direction = 0)
vol_mito = np.zeros_like(vol_dend)
vol_mito = add_shape(vol_mito, mito, [48,30,26])


# ER
er = make_cylinder(radius = 3, length = 90, direction = 0)
vol_er = np.zeros_like(vol_dend)
vol_er = add_shape(vol_er, er, [48,30,7])
vol_er = add_shape(vol_er, er, [48,30,45])


# Arrangement
vol_dend = lmpad(vol_dend)
vol_mito = lmpad(vol_mito)
vol_psd  = lmpad(vol_psd)
vol_er   = lmpad(vol_er)

vol_not_er   = np.logical_not( vol_er )
vol_not_mito = np.logical_not( vol_mito )
vol_dend_not_mito_not_er = vol_dend ^ (vol_mito | vol_er)


# Create surface
xyzpitch = 0.02
m = {}

dend = CreateSurface(vol_dend, xyzpitch, num_smoothing = 5, method_smoothing = 'laplacian')
m['bound vertices']        = dend.vertices
m['bound faces']           = dend.faces
m['bound faces in volume'] = dend.get_surface_to_volume()
m['bound faces in volume'] *= vol_dend_not_mito_not_er
face_id_psd = dend.get_face_ids_inside(vol_psd)
m['psd faces in volume']   = dend.get_surface_to_volume(face_id_psd)
m['psd faces in volume']   *= vol_dend_not_mito_not_er
m['face id psd']           = face_id_psd

mito = CreateSurface(vol_not_mito, xyzpitch)
m['mito vertices']         = mito.vertices
m['mito faces']            = mito.faces

er = CreateSurface(vol_not_er, xyzpitch)
m['er vertices']           = er.vertices
m['er faces']              = er.faces

m['unit length (um)']      = xyzpitch

tot_volume = np.zeros_like(vol_dend, dtype='uint8')
tot_volume[vol_dend > 0] = 1
tot_volume[vol_mito>0]   = 2
tot_volume[vol_er>0]     = 3
m['volume'] = tot_volume

# Save
filename = 'models/ball_and_stick.h5'
os.makedirs('models', exist_ok=True)

with h5py.File(filename,'w') as w:
	for k, v in m.items():
		w.create_dataset(k, data=v)