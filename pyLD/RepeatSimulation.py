import numpy as np
import h5py
import os, sys, errno
import configparser

class RepeatSimulation:

	def __init__(self, config_filename):
	
		if not os.path.exists(config_filename):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), config_filename)
	
		self.c = configparser.ConfigParser()
		self.c.read(config_filename, encoding='utf-8')
		pass

	def exec(self):

		for i in 

		with h5py.File(filename_prerun,'r') as f:
		    TimePoints = f['Simulations']['0000001']['Lattice'].keys()
		    TimePoints.sort()
		    print('Time ID: ', TimePoints[-1])
		    Lattice = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
		    SpeciesCount = f['Simulations']['0000001']['SpeciesCounts'][-1,:]



		with h5py.File(filename_stimrun,'a') as g:
		    g['Model']['Diffusion']['Lattice'][()] = Lattice
		    g['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount


		shutil.copy(filename_lm, filename)


