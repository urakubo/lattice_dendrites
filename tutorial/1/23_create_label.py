import sys, os
import numpy as np
import h5py
from skimage import morphology
from tut1_functions import add_shape


input_morpho_filename = 'models/ball_and_stick.h5'
output_label_filename = 'models/labels_ball_and_stick.h5'


print('Load morpho file')
with h5py.File( input_morpho_filename,'r' ) as f:
    vol_dend_not_mito_not_er = f['dendrite not mitochondrion not ER'][()]


print('Label spine')
spine_head   = morphology.ball(radius = 12)
label_volume = np.zeros_like(vol_dend_not_mito_not_er)
label_volume = add_shape(label_volume, spine_head, [48,30,76])

label_ids    = np.array([1])
label_volume = (label_volume > 0) * label_ids[0]
ref_volume   = vol_dend_not_mito_not_er


print('Save label')
with h5py.File(output_label_filename, 'a') as f:
	f['label volume'] = label_volume
	f['label ids']    = label_ids
	f['ref volume']   = ref_volume