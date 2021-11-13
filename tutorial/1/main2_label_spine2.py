import sys, os
import numpy as np
from mayavi import mlab
from mayavi.api import OffScreenEngine

main_dir = os.path.abspath(os.path.dirname(sys.argv[0]))  # Dir of main
sys.path.append(os.path.join(main_dir, '..', 'github_pages'))
from pyLD import *

input_label_file  = 'labels_ball_and_stick.h5'
output_image_file = 'labels_ball_and_stick.png'

c = LoadLabeledVolume(input_label_file)

mlab.figure(bgcolor=(1.0,1.0,1.0), size=(700,700))
mlab.view(90, 90, 300, [ 50, 30, 50 ] )

pitch = 1
for id in c.label_ids:
	vert, face,_ ,_ = create_surface(pitch, c.label_volume == id)
	vert = vert
	mlab.triangular_mesh(vert[:,0], vert[:,1], vert[:,2], face, color=tuple(np.random.rand(3))  , opacity=0.3)

vert, face, _, _ = create_surface(pitch, c.ref_volume ^ (c.label_volume > 0))
vert = vert
mlab.triangular_mesh(vert[:,0], vert[:,1], vert[:,2], face, color=(0.8,0.8,0.8)  , opacity=0.3)
mlab.savefig(output_image_file)
mlab.show()


