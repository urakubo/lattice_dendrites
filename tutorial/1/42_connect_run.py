import sys, os
import numpy as np
from pyLD import *


def event_replace(lattice, sys_param, event_param):
	i    = sys_param['i']
	time = sys_param['time']
	s    = sys_param['species']
	src  = s[event_param['src']]
	dst  = s[event_param['dst']]
	print('\nEvent at: {:g}, Current time: {:.3f}\n'.format(i, time))

	labeled_area = (sys_param['label volume'] == 1 )
	for i in range(lattice.shape[3]):
		t_lattice = lattice[:,:,:,i]
		t_lattice[labeled_area*(t_lattice == src)] = dst
	return lattice, event_param


r = ConnectRun()
r.label_volume_file = 'models/labels_ball_and_stick.h5'

r.template_lm_file  = 'models/photobleach.lm'
r.output_dir        = 'results_photobleach'
r.event_params      = {'src':'YFP', 'dst':'bleached YFP'}
r.exec_periods      = [4.0, 4.0]
r.exec_events       = [event_null, event_replace]


'''
r.template_lm_file  = 'models/Ca_dynamics.lm'
r.output_dir        = 'results_Ca_dynamics'
r.event_params      = {'src':'inactive NMDAR', 'dst':'active NMDAR'}
r.exec_periods      = [2.0, 2.0]
r.exec_events       = [event_null, event_replace]
'''

'''
GPU_ids  = '0,1'
num_GPUs = '2'
num_CPUs = '2'
r.lm_option   = ['-g', GPU_ids, '-gr', num_GPUs,'-cr', num_CPUs, '-sl','lm::rdme::MGPUMpdRdmeSolver']
'''
r.exec()
