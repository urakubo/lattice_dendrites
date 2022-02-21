import sys, os
import numpy as np
import h5py
from mayavi import mlab
from pyLD import *

input_file        = 'models/ball_and_stick.h5'
input_label_file  = 'models/labels_ball_and_stick.h5'
output_image_file = 'imgs/labels_ball_and_stick.png'


with h5py.File(input_file, 'r') as f:
	ref_volume = f['volume'][()]

with h5py.File(input_label_file, 'r') as f:
	label_volume = f['label volume'][()]
	label_ids    = f['label ids'][()]

mlab.figure(bgcolor=(1.0,1.0,1.0), size=(700,700))
mlab.view(90, 90, 300, [ 50, 30, 50 ] )

pitch = 1
for id in label_ids:
	s = CreateSurface(label_volume == id, pitch)
	color = tuple(np.random.rand(3))
	mlab.triangular_mesh(s.vertices[:,0], s.vertices[:,1], s.vertices[:,2],\
		s.faces, color=color, opacity=0.3)

d = CreateSurface( (ref_volume == 1) ^ (label_volume > 0),  pitch)
color = (0.8,0.8,0.8)
mlab.triangular_mesh(d.vertices[:,0], d.vertices[:,1], d.vertices[:,2],\
	d.faces, color=color, opacity=0.3)
mlab.savefig(output_image_file)
mlab.show()
