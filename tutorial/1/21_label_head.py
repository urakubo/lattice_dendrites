import sys, os
import numpy as np
import h5py
from skimage import morphology
from tut1_functions import add_shape
from pyLD import *

output_label_filename = 'models/labels_ball_and_stick.h5'

print('Label a spine head')
spine_head = morphology.ball(radius = 12)
vol_dend   = np.zeros((96,60,96), dtype=np.uint8)
vol_dend   = add_shape(vol_dend, spine_head, [48,30,76])
vol_dend   = lmpad(vol_dend)

print('Save label')
with h5py.File(output_label_filename, 'a') as f:
	f['label volume'] = vol_dend
	f['label ids']    = np.array([1])
