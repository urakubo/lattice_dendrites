import sys, os, errno
import numpy as np
import h5py
from pyLD import *


zpitch   = 0.04 # XX um per EM slice
xyzpitch = 0.02 # xyz-pitch (XX um) in the decoded volume

dir_dxfs = 'dxf_files'
dend     = 'contour'
er       = 'endoplasmic_reticulum'
mito     = 'mitochondrion'
psd      = 'postsynaptic_density'
save_filename = 'models/realistic_dendrite.h5'
os.makedirs('models', exist_ok=True)


print('\nCreate volumes.')
c = CreateVolumeFromReconstruct(dir_dxfs, xyzpitch, zpitch, dend)
volumes = [c.create(d) for d in [dend, er, mito, psd]]


print('\nRotate volumes on the xz plane, to obtain a minimal bounding box.')
vol_dend = volumes[0]
r = RotateVolume(vol_dend, 1)
volumes = [r.rotate(v) for v in volumes]

#print('\nSwap axes to obtain the longest z direction.)
#volumes = [v.swapaxes(1, 2) for v in volumes]


print('\nSet space as a multiple of 32 x 32 x 32 voxels.')
volumes = [lmpad(v) for v in volumes]


print('\nObtain domain volumes from the list of volumes.')
vol_dend, vol_er, vol_mito, vol_psd = volumes


print('\nSmooth edges.')
vol_dend = smooth_volume_erosion_then_dilation(vol_dend) 
vol_er   = smooth_volume_dilation_then_erosion(vol_er)
vol_psd  = smooth_volume_dilation_then_erosion(vol_psd) 


print('\nLogical operation of volumes.')
vol_not_mito = np.logical_not( vol_mito )
vol_not_er   = np.logical_not( vol_er )
vol_dend_not_mito_not_er = vol_dend ^ (vol_mito | vol_er)


print('\nCreate surface.')
dend = CreateSurface(vol_dend, xyzpitch, num_smoothing = 5, method_smoothing = 'laplacian')

m = {}
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
m['dendrite not mitochondrion not ER'] = vol_dend_not_mito_not_er

tot_volume = np.zeros_like(vol_dend, dtype='uint8')
tot_volume[vol_dend > 0]   = 1
tot_volume[vol_mito>0]     = 2
tot_volume[vol_er>0]       = 3
m['volume']                = tot_volume


with h5py.File(save_filename,'w') as w:
	for k, v in m.items():
		w.create_dataset(k, data=v)