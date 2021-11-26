import numpy as np
import h5py
import os, sys, shutil, errno
import subprocess as s
from .utils import get_species_names, get_spacing


def event_inherit(lattice, sys_param, event_param):
	"""Inherit event function for ConnectRun.
	In this event, the final state of a specified run is copied to the initial state of the next run.

	Args:
		lattice (numpy[uint8]): See "pyLD.null_event"
		sys_param (dict): See "pyLD.null_event"
		event_param (dict): Parameters for the event function that contain 

			- 'filename_prerun' (str): Filename of the prerun result.

	Returns: Tuple containing:

		- lattice: (numpy[uint8]): See "pyLD.null_event"
		- event_param: (dict): See "pyLD.null_event"
	"""
	i    = sys_param['i']
	time = sys_param['time']
	print('\nInherit event at: {:g}, Current time: {:.3f}\n'.format(i, time))

	with h5py.File(event_param['filename_prerun'],'r') as f:
		TimePoints = list( f['Simulations']['0000001']['Lattice'].keys() )
		TimePoints.sort()
		lattice_prerun = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]

	lattice = lattice_prerun
	return lattice, event_param


def event_null(lattice, sys_param, event_param):
	"""Null event function for ConnectRun. Any operation will not be executed.
	Users can modify the initial state of molelcules of each run (lattice; 4D array, 3D volume + 16 species space).

	Args:
		lattice (numpy[uint8]): Lattice space (4D array, 3D space plus 16 slots)
		sys_param (dict): System parameters that contains

			- 'i' (int): Exec id,
			- 'time' (float): Start time
			- 'species' (dict): Molecular names that have their own ids
			- 'label volume' (numpy[int]): Label volume if specified (3D array, optional)
			- 'label ids' (numpy[int]): Label ids if specified (1D array, optional)

		event_param (dict):  It contains user-defined parameters for event function

	Returns: Tuple containing:

		- lattice: (numpy[uint8]): Lattice space (4D array, 3D space plus 16 slots)
		- event_param: (dict): User-defined values can be passed to the next event
	"""

	i    = sys_param['i']
	time = sys_param['time']
	print('\nNull event at: {:g}, Current time: {:.3f}\n'.format(i, time))
	return lattice, event_param


def event_activate(lattice, sys_param, event_param):
	"""Activate event function for ConnectRun.
	In this event, a specified species of molecules ( event_param['source species'] )
	is replaced with another species ( event_param['destination species'] ).

	Args:
		lattice (numpy[uint8]): See "pyLD.null_event"
		sys_param (dict): See "pyLD.null_event"
		event_param (dict): Parameters for the event function that contain 

			- 'source species' (str): Source species
			- 'destination species' (str): Destination species

	Returns: Tuple containing:

		- lattice: (numpy[uint8]): See "pyLD.null_event"
		- event_param: (dict): See "pyLD.null_event"
	"""
	i    = sys_param['i']
	time = sys_param['time']
	print('\nActivate event at: {:g}, Current time: {:.3f}\n'.format(i, time))

	s    = sys_param['species']
	label_volume = sys_param['label volume']

	src  = event_param['source species']
	dst  = event_param['destination species']
	# prob = event_param['probability']
	# ids  = event_param['target label ids']
	# (np.random.rand(n) < prob)

	# lattice, source_molecule, dest_molecule, probability, domain = None, label_volume=):
	lattice[lattice == s[src]] = s[dst]
	return lattice, event_param


class ConnectRun:
	"""Connect multiple runs of simulation. Output of each run is set to be an inital state of the next run.
	Users can modify the initial state (4D array; 3D volume + 16 species space through user-defined event functions.
	All the argments below are assigned as instance variables of this class.
	Thus, users can set/modify them anytime before "exec".

	Args:
		template_lm_file (str): Template lm file
		lm_option (list): Optional arguments passing to lm
		exec_periods (list[float]): Simulation period in each run
		exec_events (list[obj]): Event function to modify the lattice space before each run
		event_params (dict/list[dict]/tuple[dict]): Parameters passing to event function
		output_dir (str): Directory that stores simulation results
		output_prefix (str): Prefix of output lm filenames that stores simulation results
		output_num_zero_padding (int): Number of zero padding in output lm filenames
		start_id (int): Start id for the output filename
		label_volume_file (str): labeled volume file (optional)

	Returns:
		(pyLD.ConnectRun): ConnectRun object that contains the above instance variables
	"""
	def __init__(self,
		template_lm_file = None,
		lm_option     = ['-r', '1', '-sp', '-sl','lm::rdme::MpdRdmeSolver'],
		exec_periods  = [1],
		exec_events   = [event_null],
		event_params  = None,
		output_dir    = 'results',
		output_prefix = '',
		output_num_zero_padding = 4,
		label_volume_file = None,
		start_id = 0
		):

		self.template_lm_file = template_lm_file
		self.lm_option        = lm_option
		self.exec_periods     = exec_periods
		self.exec_events      = exec_events
		self.event_params     = event_params
		self.output_dir       = output_dir
		self.output_prefix    = output_prefix
		self.output_num_zero_padding = output_num_zero_padding
		self.label_volume_file= label_volume_file
		self.start_id         = start_id


	def exec(self):
		"""Execute repeat runs,

		Args:

		Returns: bool 
			(bool): True if succeeded. Also simulation results are stored in lm files in output_dir.
		"""
		self.exec_periods     = list(self.exec_periods)
		self.exec_events      = list(self.exec_events)
		if self.template_lm_file==None or not os.path.isfile(self.template_lm_file):
			raise ValueError('No template lm file.')
		elif len(self.exec_periods) != len(self.exec_events):
			raise ValueError('Num of exec_periods must be the same as the num of exec_events.')
		elif not isinstance(self.event_params, dict) and (len(self.event_params) != len(self.exec_events)):
			raise ValueError('usr_params must be dict or a list/tuple of dict that has the same length with exe_events')

		# Set system params
		sys_param = {}
		sys_param['time'] = 0.0
		sys_param['species'] = get_species_names(self.template_lm_file)
		if self.label_volume_file != None:
			sys_param = self._load_label_file(sys_param)

		os.makedirs(self.output_dir, exist_ok=True)
		filename = ''
		filename_prerun = ''
		sys_param['i'] = self.start_id
		for i, (period, event) in enumerate(zip(self.exec_periods, self.exec_events)):
			# Copy results from a previous run or inits from the orignal template.
			if i > 0:
				sys_param['time'] = sys_param['time'] + self.exec_periods[i-1]
				with h5py.File(filename_prerun,'r') as f:
				    TimePoints = list( f['Simulations']['0000001']['Lattice'].keys() )
				    TimePoints.sort()
				    lattice = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
				    species_count = f['Simulations']['0000001']['SpeciesCounts'][-1,:]
			else:
				with h5py.File(self.template_lm_file,'r') as f:
					lattice = f['Model']['Diffusion']['Lattice'][()]
					species_count = f['Model']['Reaction']['InitialSpeciesCounts'][()]

			# Execute an event to change the lattice
			if isinstance(self.event_params, dict):
				lattice, self.event_params = event(lattice, sys_param, self.event_params )
			else:
				lattice, _ = event(lattice, sys_param, self.event_params[i] )
			uniq, count = np.unique(lattice, return_counts=True)
			u_ids = np.where(uniq > 0)
			uniq  = uniq[u_ids]
			count = count[u_ids]
			# print('uniq : ', uniq)
			# print('count: ', count)
			species_count         = np.zeros_like(species_count)
			species_count[uniq-1] = count

			# Copy a lm file from the original template
			filename = self.output_prefix + str(sys_param['i']).zfill(self.output_num_zero_padding)+'.lm'
			filename = os.path.join(self.output_dir, filename)
			shutil.copy(self.template_lm_file, filename)

            # Modify the lm file
			with h5py.File(filename,'a') as f:
				period_s = str(period).ljust(4)[:4]
				f['Parameters'].attrs['maxTime'] = np.string_(period_s)
				f['Model']['Diffusion']['Lattice'][()] = lattice
				f['Model']['Reaction']['InitialSpeciesCounts'][()] = species_count

            # Execute the modified lm file
			command = ['lm', '-f', filename]
			command.extend(self.lm_option)
			print(' '.join(command))
			s.call(command)
			filename_prerun = filename
			sys_param['i']  += 1

		return True


	def _load_label_file(self, sys_param):
		with h5py.File(self.label_volume_file, 'r') as f:
			sys_param['label volume'] = f['label volume'][()]
			sys_param['label ids']    = f['label ids'][()]
		return sys_param

