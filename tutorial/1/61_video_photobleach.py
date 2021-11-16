import os
import numpy as np
from mayavi import mlab
from mayavi.api import OffScreenEngine
from pyLD import *
from tut1_functions import show_dendrite
import subprocess as s

# mlab.options.offscreen = True
morpho_file = "models/ball_and_stick.h5"
lm_files    = ['results_photobleach/0000.lm'  ,'results_photobleach/0001.lm']
output_dir  = "imgs_for_video"
os.makedirs(output_dir, exist_ok=True)

# Show spiny dendrite
mlab = show_dendrite(morpho_file)


def event(lattice, sys_param, event_param):
	id   = sys_param['id']
	time = sys_param['time']
	s    = sys_param['species']
	print('Event at: {:g}, time: {:.3f}'.format(id, time))

	YFP   = []
	for i in range(lattice.shape[3]):
		YFP.extend( np.flatnonzero(lattice[:,:,:,i] == s['YFP']).tolist() )
	YFP = np.unravel_index(YFP, lattice[:,:,:,0].shape )
	plot_points = mlab.points3d(YFP[0], YFP[1], YFP[2], color=(0.8,0.8,0), scale_factor=1.5,line_width=0.1)
	text = mlab.text(0.2,0.05,"{0:.3f} s".format(time) ,color=(0,0,0), width=0.3)

	mlab.view( 90, 90, 300, [ 50, 30, 50 ] )
	fname = os.path.join(output_dir, str(id).zfill(4)+ '.png')
	mlab.savefig(fname)
	plot_points.remove()
	text.remove()

	return event_param

print('\nGenerate images\n')
c = ConnectAnalysis(lm_files)
c.start_time  = -4
c.event       = event
c.event_param = {}
c.exec()


print('\nConnect them into a video\n')
ffname = os.path.join(output_dir, '%04d'+ '.png')
com = ['ffmpeg','-r', '10', '-i', ffname,'-pix_fmt', 'yuv420p', 'photobleach.mp4']
print(' '.join(com))
s.call(com)