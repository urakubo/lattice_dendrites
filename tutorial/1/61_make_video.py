import os
import shutil
import numpy as np
from mayavi import mlab
import subprocess as s
from pyLD import *
from tut1_functions import *

# FRAP
targ_name   = 'YFP'
color       = (1,1,0)
lm_files    = ['results_photobleach/0000.lm'  ,'results_photobleach/0001.lm']
start_time  = -4

'''
# Ca2+ influx via NMDARs
targ_name   = 'Ca'
color       = (0,0,1)
lm_files    = ['results_Ca_dynamics/0000.lm'  ,'results_Ca_dynamics/0001.lm']
start_time  = -2
'''

morpho_file = "models/ball_and_stick.h5"
output_dir  = "imgs_for_video"
if os.path.isdir(output_dir):
	shutil.rmtree(output_dir)
os.makedirs(output_dir, exist_ok=True)


# Show spiny dendrite
mlab = show_dendrite(morpho_file)


# Event to make an image
def event(lattice, sys_param, event_param):
	id   = sys_param['id']
	time = sys_param['time']
	s    = sys_param['species']
	print('Event at: {:g}, time: {:.3f}'.format(id, time))
	m    = get_loc_molecule(lattice, s[targ_name])
	plot_points = mlab.points3d(m[0], m[1], m[2],\
		color = color,\
		scale_factor = 1.5)
	text = mlab.text(0.2,0.05,"{0:.3f} s".format(time), color=(0,0,0), width=0.3)

	mlab.view( 90, 90, 300, [ 50, 30, 50 ] )
	fname = os.path.join(output_dir, str(id).zfill(4)+ '.png')
	mlab.savefig(fname)
	plot_points.remove()
	text.remove()

	return event_param


print('\nGenerate images')
c = ConnectAnalysis(lm_files)
c.start_time  = start_time
c.event       = event
c.event_param = {}
c.exec()


print('\nConnect images into a video')
ffname = os.path.join(output_dir, '%04d.png')
com = ['ffmpeg','-r', '10', '-i', ffname,\
       '-vf','scale=trunc(iw/2)*2:trunc(ih/2)*2',\
       '-pix_fmt', 'yuv420p', targ_name+'.mp4']
print(' '.join(com))
s.call(com)