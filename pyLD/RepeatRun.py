import numpy as np
import h5py
import os, sys, shutil, errno
import subprocess as s
import configparser
import json

ABSOLUTE_PATH    = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INI_FILE = os.path.join(ABSOLUTE_PATH, 'repeat_default.ini')

class RepeatRun:

	def __init__(self, config_filename):
		if not os.path.exists(config_filename):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), config_filename)

		# Load
		# https://stackoverflow.com/questions/6107149/how-can-i-set-default-values-for-safeconfigparser
		# self.c = configparser.ConfigParser()
		c = configparser.SafeConfigParser()
		c.read(DEFAULT_INI_FILE, encoding='utf-8')
		c.read(config_filename , encoding='utf-8')

		self.template_lm_file = c['Input']['template_lm_file']
		self.directory        = c['Output']['directory']
		self.lmfile_prefix    = c['Output']['lmfile_prefix']
		self.num_zero_padding = c.getint('Output', 'num_zero_padding')
		self.num_repeats      = c.getint('Fixed repeat', 'number')
		self.period           = c.getfloat('Fixed repeat', 'period')
		self.periods          = np.array( json.loads(c.get("Variable repeat","periods")) )
		self.command          = json.loads(c.get("Simulation setting","command").replace("'", '"'))

	def exec(self):

		if not os.path.isfile(self.template_lm_file):
			raise ValueError('No lm file.')
		elif (np.all(self.periods > 0)):
			periods = self.periods
		elif (self.period > 0) and (self.num_repeats > 0):
			periods = np.ones(self.num_repeats) * self.period
		else :
			raise ValueError('Invalid num_repeats/period/periods.')

		os.makedirs(self.directory, exist_ok=True)

		filename = ''
		filename_prerun = ''
		for i, period in enumerate(list(periods)):

			# Make lm file from a original template
			filename = os.path.join(self.directory, self.lmfile_prefix + str(i).zfill(self.num_zero_padding)+'.lm')
			shutil.copy(self.template_lm_file, filename)

			# Copy results from the previous run
			if i != 0:
				with h5py.File(filename_prerun,'r') as f:
				    TimePoints = f['Simulations']['0000001']['Lattice'].keys()
				    TimePoints.sort()
				    print('Time ID: ', TimePoints[-1])
				    Lattice = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
				    SpeciesCount = f['Simulations']['0000001']['SpeciesCounts'][-1,:]
				with h5py.File(filename,'a') as g:
					g['Model']['Diffusion']['Lattice'][()] = Lattice
					g['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount

			# Modify the lm file
			t = str(period) + '     '
			t = t[:4]
			with h5py.File(filename,'a') as f:
				f['Parameters'].attrs['maxTime'] = t


#        Lattice[Lattice == S['NR']] = S['NR_Glu']
#        SpeciesCount[S['NR_Glu']-1] = SpeciesCount[S['NR_Glu']-1]\
#                                      + SpeciesCount[S['NR']-1]
#        SpeciesCount[S['NR']-1]     = 0


# (Remove, add, replace) in (domains, surfaces)
# Manual operation

# f = h5py.File(filename_prerun[-1],'r')
# mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
# S = {}
# for i in range(len(mnames)):
#    S[mnames[i]] = i+1
# test = f['Parameters'].attrs['maxTime']
# f.close()



			command = self.command + ['-f', filename]
			print(' '.join(command))
			s.call(command)

			filename_prerun = filename

