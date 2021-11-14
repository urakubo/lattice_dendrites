import sys, os, glob, pickle
import sys, os
from pyLD import *



def activate_NMDAR(lattice, sys_param, usr_params):
	tmp = sys_param['lattice_volume']
	return lattice, usr_params


#GPU_ids  = '0,1'
#num_GPUs = '2'
#num_CPUs = '2'
#option   = ['-g', GPU_ids, '-gr', num_GPUs,'-cr', num_CPUs, '-sl','lm::rdme::MGPUMpdRdmeSolver']


exec_periods = [1.0, 1.0]
exec_events  = [activate_NMDAR, activate_NMDAR]
r = RepeatRun()


