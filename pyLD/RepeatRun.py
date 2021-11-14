import numpy as np
import h5py
import os, sys, shutil, errno
import subprocess as s
from .utils import get_species_names, get_spacing


def null_event(lattice, sys_param, usr_params):
	return lattice, usr_params


def activate(lattice, sys_param, usr_param):
	s    = sys_param['species']
	label_volume = sys_param['label volume']
	src  = usr_param['source species']
	dst  = usr_param['destination species']
	prob = usr_param['probability']
	ids  = usr_param['target label ids']

	# (np.random.rand(n) < prob)

	# lattice, source_molecule, dest_molecule, probability, domain = None, label_volume=):
	lattice[lattice == s[src]] = s[dst]
	return lattice, usr_params


class RepeatRun:
	"""Repeat simulation runs. Output of each run is set to be a inital state of the next run.
	Users can modify the initial state (4D array; 3D volume + 16 species space),
	which can be considered as event function.

	Returns:
		(pyLD.RepeatRun): RepeatRun object that contains the following instance variables:

		- 'template_lm_file' (str): Template lm file
		- 'lm_option' (list): Optional arguments passing to lm
		- 'exec_periods' (list[float]): Simulation period in each run
		- 'exec_events' (list[obj]): Event function to modify the lattice space before each run
		- 'usr_params' (dict/list[dict]): Usr parameters passing to each event function
		- 'output_dir' (str): Directory that stores simulation results
		- 'output_prefix' (str): Prefix of output lm filenames that stores simulation results
		- 'output_num_zero_padding' (int): Number of zero padding in output lm filenames
		- 'labeled_volume_file' (str): labeled volume file (optional)
	"""
	def __init__(self,
		template_lm_file=None,
		lm_option    = ['-r', '1', '-sp', '-sl','lm::rdme::MpdRdmeSolver'],
		exec_periods = [1],
		exec_events  = [null_event],
		usr_params   = None,
		output_dir='results',
		output_prefix=None,
		output_num_zero_padding=4,
		labeled_volume_file=None
		):

		self.template_lm_file = template_lm_file
		self.lm_option        = lm_option
		self.exec_periods     = exec_periods
		self.exec_events      = exec_events
		self.usr_params       = usr_params
		self.output_dir       = output_dir
		self.output_prefix    = output_prefix
		self.output_num_zero_padding = output_num_zero_padding
		self.labeled_volume_file = labeled_volume_file


	def exec(self):
		"""Execute repeat runs,

		Args:

		Returns: bool 
			(bool): True if succeeded. Also simulation results are stored in lm files in output_dir.
		"""

		self.exec_periods     = list(self.exec_periods)
		self.exec_events      = list(self.exec_events)
		if not os.path.isfile(self.template_lm_file):
			raise ValueError('No template lm file.')
		elif len(self.exec_periods) != len(self.exec_events):
			raise ValueError('Num of exec_periods must be the same as the num of exec_events.')
		elif not isinstance(self.usr_params, dict) or (len(self.usr_params) != len(self.exec_events)):
			raise ValueError('usr_params must be dict or a list/tuple of dict that has the same length with exec_events.')

		# Set system params
		sys_param = {}
		sys_param['time'] = 0.0
		sys_param['s']    = get_species_names(self.template_lm_file)
		if self.labeled_volume_file != None:
			sys_param = self._load_label_file(sys_param)

		os.makedirs(self.output_dir, exist_ok=True)
		filename = ''
		filename_prerun = ''
		for i, (period, event) in enumerate(zip(self.exec_periods, self.exec_events)):

			# Make lm file from a original template
			file = self.output_prefix + str(i).zfill(self.output_num_zero_padding)+'.lm'
			filename = os.path.join(self.output_dir, file)
			shutil.copy(self.template_lm_file, filename)

			# Copy results from the previous run
			sys_param['i'] = i
			if i > 0:
				sys_param['time'] += self.exec_periods[i-1]
				with h5py.File(filename_prerun,'r') as f:
				    TimePoints = f['Simulations']['0000001']['Lattice'].keys()
				    TimePoints.sort()
				    lattice = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
				    species_count = f['Simulations']['0000001']['SpeciesCounts'][-1,:]
			else:
				with h5py.File(filename,'r') as f:
					lattice = f['Model']['Diffusion']['Lattice'][()]
					species_count = f['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount

			# Modify the lm file
			if isinstance(self.usr_params, dict):
				lattice, self.usr_params = event(lattice, sys_param, self.usr_params )
			else:
				lattice, _ = event(lattice, sys_param, self.usr_params[i] )

			with h5py.File(filename,'a') as f:
				f['Parameters'].attrs['maxTime'] = str(period).ljust(4)[:4]
				f['Model']['Diffusion']['Lattice'][()] = lattice
				f['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount

			command = ['lm', '-f', filename]
			command.extend(self.lm_option)
			print(' '.join(command))
			s.call(command)
			filename_prerun = filename
		return True

	def _load_label_file(self, sys_param):
		with h5py.File(self.labeled_volume_file, 'r') as f:
			sys_param['label volume'] = f['label volume'][()]
			sys_param['label ids']    = f['label ids'][()]
		return sys_param


