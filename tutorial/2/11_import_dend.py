import sys, os, errno
import numpy as np
import h5py
from pyLD import *


zpitch   = 0.04 # XX um per EM slice
xyzpitch = 0.02 # XX um per xyz-voxel in the volume

folder = r'dxf'
dend   = 'dend_cyan01'
er     = 'dend_cyan01_ER'
mito   = 'dend_cyan01_Mito'
psd    = 'dend_cyan01_PSD'
filename     = 'dend01.h5'
annot_folder = 'annot_dend01'


# Create volume
c = CreateVolumeFromReconstruct(folder, xyzpitch, zpitch, dend)
volumes = [c.create(d) for d in [dend, er, mito, psd]]


# Rotate volumes into a minimal bounding box of the xz plane
vol_dend = volumes[0]
r = RotateVolume(vol_dend, 1)
volumes = [r.rotate(v) for v in volumes]


# Swap axes for the longest z direction, depending on dendrites.
#volumes = [v.swapaxes(1, 2) for v in volumes]


# Pad spaces for lm simulation
volumes = [lmpad(v) for v in volumes]


# Extract domains from a list of volumes
vol_dend, vol_er, vol_mito, vol_psd = volumes


# Smooth volumes
vol_dend = smooth_volume_erosion_then_dilation(vol_dend) 
vol_er   = smooth_volume_dilation_then_erosion(vol_er)
vol_psd  = smooth_volume_dilation_then_erosion(vol_psd) 


# Logical operation of volumes
vol_not_mito = np.logical_not( vol_mito )
vol_not_er   = np.logical_not( vol_er )
vol_dend_not_mito_not_ER = vol_dend ^ (vol_mito | vol_er)


# Create surface
bound_verts, bound_faces, bound_area_per_face, bound_areas, id_face_psd = \
	create_surface(xyzpitch, vol_dend, PSD = vol_psd, num_smoothing = 5, method_smoothing = 'laplacian')
mito_verts, mito_faces, _, mito_areas = create_surface(xyzpitch, vol_not_mito)
er_verts, er_faces, _, er_areas = create_surface(xyzpitch, vol_not_er)

mito_areas  = mito_areas  * vol_dend_not_mito_not_er
er_areas    = er_areas    * vol_dend_not_mito_not_er
bound_areas = bound_areas * vol_dend_not_mito_not_er


# Save
with h5py.File(filename,'w') as w:
    w['unit length per voxel (um)'] = xyzpitch
    w['dendrite']                   = vol_dend.astype(np.uint8)
    w['PSD']                        = vol_psd.astype(np.uint8)
    w['mitochondrion']              = vol_mito.astype(np.uint8)
    w['er']                         = vol_er.astype(np.uint8)
    w['dendrite not mitochondrion not ER'] = vol_dend_not_mito_not_er

    w['boundary areas in volume']   = bound_areas
    w['boundary vertices']          = bound_verts
    w['boundary faces']             = bound_faces
    w['PSD ids in boundary faces']  = id_face_psd

    w['mitochondrion areas in volume'] = mito_areas
    w['mitochondrion vertices']        = mito_verts
    w['mitochondrion faces']           = mito_faces

    w['er areas in volume'] = er_areas
    w['er vertices']        = er_verts
    w['er faces']           = er_faces


# Save UNI-EM annot
bound_color = [192,192,192]
mito_color  = [255,255,152]
er_color    = [179,255,179]
surfaces = {1: [bound_verts, bound_faces, bound_color],\
			2: [mito_verts, mito_faces, mito_color],\
			3: [er_verts, er_faces, er_color]}
save_uniem_annotator(annot_folder, xyzpitch, (vol_dend+vol_mito+vol_er*2).astype('uint16'), surfaces)