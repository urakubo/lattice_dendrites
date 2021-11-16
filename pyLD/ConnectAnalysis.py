import numpy as np
import h5py
import os, sys
from .utils import get_species_names, get_spacing, num_to_uM


class ConnectAnalysis:
	"""Enable user-defined analyses on the time developments of molecules from a series of lm files.
	
	Args:
		lm_files (str / list[str] / tuple[str]): Simulated lm files.

	Returns:
		(pyLD.ConnectAnalysis): ConnectAnalysis object
	"""

	def  __init__(self, lm_files):
		# Check arguments
		if isinstance(lm_files, str):
		    lm_files = [lm_files]
		elif isinstance(lm_files, list) | isinstance(lm_files, tuple) :
			pass
		else:
			raise ValueError('lm_files must be str, list, or tuple.')
		self.lm_files = lm_files

		# Timepoints
		for i, fname in enumerate(self.lm_files):
		    with h5py.File(fname,'r') as f:
		        t = f['Simulations']['0000001']['LatticeTimes'][()]
		        n = f['Simulations']['0000001']['SpeciesCounts'][()]
		    if i == 0:
		        self.timepoints = t
		        self.numbers    = n
		    else:
		        self.timepoints = np.hstack( (self.timepoints, t[1:]+self.timepoints[-1]) )
		        self.numbers    = np.concatenate([self.numbers , n[1:,:] ], axis=0)

		self.species  = get_species_names(self.lm_files[0])
		self.label_volume_file = None
		self.event             = None
		self.event_param       = None
		self.start_time        = 0.0

	def exec(self):
		"""Execute repeats

		Args:

		Returns: bool 
			(bool): True if succeeded. Also simulation results are stored in lm files in output_dir.
		"""
		sys_param = {}
		sys_param['species']    = self.species
		sys_param['timepoints'] = self.timepoints
		sys_param['numbers']    = self.numbers
		if self.label_volume_file != None:
			sys_param = self._load_label_file(sys_param)

		id = 0
		for lm_file in self.lm_files:
			self.start_time, id = self._process_each_lm_file(lm_file, sys_param, self.start_time, id)


	def _process_each_lm_file(self, lm_file, sys_param, start_time, id):

		# Time frames
		print('file :', lm_file)
		with h5py.File(lm_file,'r') as file:
			timepoints = list( file['Simulations']['0000001']['LatticeTimes'][()] )
			frames = [key for key in file['Simulations']['0000001']['Lattice'].keys()]
			frames.sort()

			if id != 0:
				timepoints = timepoints[1:]
				frames     = frames[1:]

			for t, f in zip(timepoints, frames):
				lattice = file['Simulations']['0000001']['Lattice'][f][:,:,:,:]
				sys_param['time'] = start_time + t
				sys_param['id']   = id
				self.event(lattice, sys_param, self.event_param)
				id += 1
		start_time = start_time + timepoints[-1]
		return start_time, id


	def _load_label_file(self, sys_param):
		with h5py.File(self.label_volume_file, 'r') as f:
			sys_param['label volume'] = f['label volume'][()]
			sys_param['label ids']    = f['label ids'][()]
		return sys_param

