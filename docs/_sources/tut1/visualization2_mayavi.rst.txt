=======================
Visualization 2: mayavi
=======================

まずは、好きなディレクトリに移動して、sphinx-quickstartを実行します。

以下では report というドキュメント名を指定して実行しています。


.. code-block:: python
	:linenos:

	# Generate a spiny dendrite.

	import sys, os
	import numpy as np
	import h5py
	from skimage import morphology
	from pyLD import *

	def make_cylinder(radius, length, direction = 2):
		arrays = [morphology.disk(radius) for _ in range(length)]
		cylinder = np.stack(arrays, axis=0)
		cylinder = cylinder.swapaxes(0, direction)
		return cylinder

	def add_shape(volume, object, loc_center):
		s = np.array(object.shape)
		c = np.floor(s/2).astype(int)
		b = loc_center - c
		e = b + s
		volume[b[0]:e[0], b[1]:e[1], b[2]:e[2] ] += object
		volume = (volume > 0).astype(np.uint8)
		return volume


	# Spiny dendrite
	spine_head = morphology.ball(radius = 12)
	spine_neck = make_cylinder(radius = 5, length = 50, direction = 2)
	dendrite   = make_cylinder(radius = 25, length = 90, direction = 0)


	vol_dend = np.zeros((100,60,110), dtype=np.uint8)
	vol_dend = add_shape(vol_dend, spine_head, [50,30,105-25])
	vol_dend = add_shape(vol_dend, spine_neck, [50,30,105-50])
	vol_dend = add_shape(vol_dend, dendrite  , [50,30,105-75])


	# PSD
	psd = morphology.ball(radius = 7)
	vol_psd = add_shape(np.zeros_like(vol_dend), psd, [50,30,105-12])


	# Mito
	mito = make_cylinder(radius = 10, length = 90, direction = 0)
	vol_mito = np.zeros_like(vol_dend)
	vol_mito = add_shape(vol_mito, mito, [50,30,105-75])


	# ER
	er = make_cylinder(radius = 3, length = 90, direction = 0)
	vol_er = np.zeros_like(vol_dend)
	vol_er = add_shape(vol_er, er, [50,30,105-61])
	vol_er = add_shape(vol_er, er, [50,30,105-89])


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

	bound_verts, bound_faces, bound_area_per_face, bound_areas, id_face_psd = \
		create_surface(xyzpitch, vol_dend, PSD = vol_psd, num_smoothing = 5, method_smoothing = 'laplacian')

	mito_verts, mito_faces, mito_area_per_face, mito_areas = create_surface(xyzpitch, vol_not_mito)
	er_verts, er_faces, er_area_per_face, er_areas         = create_surface(xyzpitch, vol_not_er)

	mito_areas  = mito_areas  * vol_dend_not_mito_not_er
	er_areas    = er_areas    * vol_dend_not_mito_not_er
	bound_areas = bound_areas * vol_dend_not_mito_not_er

	# Save UNI-EM annot
	annot_folder = 'annot_ball_and_stick'
	bound_color = [192,192,192]
	mito_color  = [255,255,152]
	er_color    = [179,255,179]
	surfaces = {1: [bound_verts, bound_faces, bound_color],\
				2: [mito_verts, mito_faces, mito_color],\
				3: [er_verts, er_faces, er_color]}
	save_uniem_annotator(annot_folder, xyzpitch, (vol_dend+vol_mito+vol_er*2).astype('uint16'), surfaces)


	# Save

	filename = 'ball_and_stick.h5'

	with h5py.File(filename,'w') as w:
	    w['unit length per voxel (um)'] = xyzpitch
	    w['dendrite'] 					= vol_dend.astype(np.uint8)
	    w['PSD']      					= vol_psd.astype(np.uint8)
	    w['mitochondrion']      		= vol_mito.astype(np.uint8)
	    w['er']	      					= vol_er.astype(np.uint8)
	    w['dendrite not mitochondrion not ER']  = vol_dend_not_mito_not_er

	    w['boundary areas in volume']   = bound_areas
	    w['boundary vertices']      	= bound_verts
	    w['boundary faces']        		= bound_faces
	    w['PSD ids in boundary faces'] 	= id_face_psd

	    w['mitochondrion areas in volume'] = mito_areas
	    w['mitochondrion vertices']      = mito_verts
	    w['mitochondrion faces']         = mito_faces

	    w['er areas in volume'] = er_areas
	    w['er vertices']        = er_verts
	    w['er faces']           = er_faces


今回は、index.rstの中身は必要ないので、一旦全て削除します。
