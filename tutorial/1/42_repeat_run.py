import sys, os
import numpy as np
from pyLD import *


def remove(lattice, sys_param, event_param):
	i    = sys_param['i']
	time = sys_param['time']
	print('\nRemove event at: {:g}, Current time: {:.3f}\n'.format(i, time))

	target_volume = (sys_param['label volume'] == event_param['target label id'] )
	tmp_lattice = lattice
	for i in range(lattice.shape[3]):
		tmp_lattice = lattice[:,:,:,i]
		tmp_lattice[target_volume] = 0
		lattice[:,:,:,i] = tmp_lattice
	return lattice, event_param


r = RepeatRun()
r.exec_periods = [4.0, 4.0]
r.exec_events  = [null_event, remove]
r.event_params = {'species': 'YFP', 'target label id': 1}
r.template_lm_file    = 'models/ball_and_stick_photobleach.lm'
r.labeled_volume_file = 'models/labels_ball_and_stick.h5'
r.output_dir          = 'results_photobleach'

#GPU_ids  = '0,1'
#num_GPUs = '2'
#num_CPUs = '2'
#r.lm_option   = ['-g', GPU_ids, '-gr', num_GPUs,'-cr', num_CPUs, '-sl','lm::rdme::MGPUMpdRdmeSolver']
r.exec()