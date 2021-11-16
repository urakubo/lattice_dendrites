import numpy as np
import h5py
import os, sys
from .utils import get_species_names, get_spacing, num_to_uM

NUMPY_INTEGERS = [ np.int8, np.int16, np.int32, np.int64, \
	np.uint8, np.uint16, np.uint32, np.uint64]


class GetLabeledConcs:
	"""Get time series of moleuclar numbers/concentrations within labeled volumes from LM simulation result.

	Args:
		lm_file (str): Target lm file
		monitor_species (str): Target molecule to monitor

	Returns:
		(pyLD.GetLabeledConcs): GetLabeledConcs object
	"""
	def  __init__(self, lm_file = None, monitor_species = None):

		self.lm_file = lm_file
		self.monitor_species = monitor_species

		self.timepoints = []
		self.concs = []
		self.numbers = []
		self.label_ids = []
		self.label_volume = None


	def _count_molecules(self, particles, targ_spine_labels, s):
		num_molecules_spine = {}
		for targ in s.keys():
			num_molecules_spine[targ] = [0 for i in range(len(targ_spine_labels))]

		for j in range(particles.shape[3]):
			p = np.ravel(particles[:,:,:,j])
			for i, targ_spine_label in enumerate(targ_spine_labels):
			    pp = p[targ_spine_label]
			    for targ in s.keys():
			        num_molecules_spine[targ][i] += np.count_nonzero( pp == s[targ] ) 
		return num_molecules_spine


	def save(self, filename):
		""" Save time series of moleuclar numbers/concentrations within labeled volumes.
		Args:
			filename (none/str): hdf5 filename.

		Returns:
			(file): Following variables are saved into a hdf5 file.:

			- timepoints (numpy[float]): Timepoints (s)
			- species (numpy[obj]): Names of species
			- ids_label (numpy[int]): Labels in the volume. The numpy array has [label1, label2, ...].
			- number (numpy[int]): Time series of the numbers of molecules in labels (3D array; [Timepoints, species, label_ids]).
			- conc in uM (numpy[float]): Time series of the concentrations of molecules in labels (unit: uM) (3D array; [Timepoints, species, label_ids]).
		"""
		if not isinstance(filename, str):
			raise ValueError('filename is not str.')
		elif (self.timepoints == []) or (self.concs == []) or (self.numbers == []) or (self.label_ids == []):
			raise ValueError('No labeled volume info.')

		data_numbers = []
		data_concs  = []
		species = []
		for targ in self.s.keys():
			species.append(targ)
			data_numbers.append(self.numbers[targ])
			data_concs.append(self.concs[targ])

		species      = np.array(species, dtype=h5py.special_dtype(vlen=str))
		data_numbers = np.swapaxes(np.array(data_numbers),0,1)
		data_concs   = np.swapaxes(np.array(data_concs),0,1)

		print('Savefile : ', filename)
		with h5py.File(filename, 'w') as f:
			f['timepoints'] = self.timepoints
			f['label_ids']  = self.label_ids
			f['species']    = species
			f['numbers']    = data_numbers
			f['concs in uM'] = data_concs


	def load_label_volume(self, label_file):

		"""Load a label volume from a label volume file.

		Args:
			label_file (str): Filename of the label volume file.

		Returns:
			(bool): True or False
		"""
		if isinstance(label_file, str) :
			with h5py.File(label_file,'r') as f:
				self.label_volume = f['label volume'][()]
		else:
			raise ValueError('label_file must be str.')

		return True


	def exec(self):
		"""Get time series of moleuclar numbers/concentrations within labeled volumes from LM simulation result.
		label_volume must be specified, which also can be obtained by "load_label_volume".

		Returns:
			(pyLD.GetLabeledConcs): GetLabeledConcs object that contains the following instance variables

				- timepoints (numpy[float]): Timepoints (s)
				- concs (dict): Time series of molecular concentrations (number per labeled volume, in the unit of uM). The dict container has {'species1': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], 'Species2': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], ... }.
				- numbers (dict): Time series of the numbers of molecules of the specified molecular species. The dict container has {'species1': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], 'Species2': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], ... }.
				- label_ids (numpy[int]): Labels in the volume. The numpy array has [label1, label2, ...].
		"""

		# Check arguments
		if not isinstance(self.lm_file, str) :
			raise ValueError("lm_file must be str.")
		elif not isinstance(self.label_volume, np.ndarray) or (self.label_volume.ndim != 3) or (self.label_volume.dtype not in NUMPY_INTEGERS):
			raise ValueError('label_volume must be a integer 3D np.ndarray.')

		self.s = get_species_names(self.lm_file)
		self.spacing = get_spacing(self.lm_file)

		label_ids, nums_label_voxels = np.unique(self.label_volume, return_counts=True)
		label_ids          = label_ids[1:] #  0 was removed.
		nums_label_voxels  = nums_label_voxels[1:]

		targ_spine_labels = []
		labels_flat = np.ravel(self.label_volume)
		for label_id in label_ids:
		    targ_spine_label = np.where( (labels_flat == label_id) )
		    targ_spine_labels.append( targ_spine_label )

		# Time frames
		with h5py.File(self.lm_file,'r') as file:
			timepoints = file['Simulations']['0000001']['LatticeTimes'][()]
			frames = [key for key in file['Simulations']['0000001']['Lattice'].keys()]
		frames.sort()

		# Monitor messages
		if self.monitor_species != None:
		    print('file :', self.lm_file)
		    print('Label ids : ', self.label_ids)
		    print('Species    : ', ', '.join(list(self.s.keys())))


		# Get molecules in labeled volumes at each timepoint
		numbers = {}
		for Targ in self.s.keys():
		    numbers[Targ] = []

		for f in frames:
			#
			num_molecules_time_i = {}
			for Targ in self.s.keys():
			    num_molecules_time_i[Targ] = [0 for i in range(max(label_ids)+1)]

			#
			with h5py.File(self.lm_file,'r') as file:
			    particles = file['Simulations']['0000001']['Lattice'][f][:,:,:,:]

			#
			num_molecules_time_i = self._count_molecules(particles, targ_spine_labels, self.s)
			if isinstance(self.monitor_species, str) and self.monitor_species in self.s.keys():
			    print('Monitor species:', self.monitor_species,', Number_at_time_i:', num_molecules_time_i[self.monitor_species] )

			#
			for Targ in self.s.keys():
				numbers[Targ].append(num_molecules_time_i[Targ])

		concs = {}
		for Targ in self.s.keys():
			numbers[Targ] = np.array(numbers[Targ])
			concs[Targ] = num_to_uM(numbers[Targ], nums_label_voxels, self.spacing)

		self.timepoints = timepoints
		self.concs = concs
		self.numbers = numbers
		self.label_ids = label_ids


