import numpy as np
import h5py
import os, sys
from .utils import get_species_names, get_spacing, num_to_uM

NUMPY_INTEGERS = [ np.int8, np.int16, np.int32, np.int64, \
	np.uint8, np.uint16, np.uint32, np.uint64]


class GetLabeledConcs:
	"""Get time series of moleuclar numbers/concentrations within labeled volumes from LM simulation result.

	Args:
		sim (obj): RDMESimulation object
		volume (numpy[int]): 3D array that specifies volume_ids
		domains (dict): {'domain name', volume_id}
		surfaces (dict): {'surface_name', numpy[float]} : Surface areas in voxel space (3D array)

	Returns:
		(pyLD.GetLabeledConcs): GetLabeledConcs object
	"""
	def  __init__(self,lm_filename, monitor_species = None):

		# Check arguments
		if not isinstance(lm_filename, str) :
			raise ValueError("lm_filename must be str.")

		self.lm_filename = lm_filename
		self.monitor_species = monitor_species
		self.s = get_species_names(lm_filename)
		self.spacing = get_spacing(lm_filename)

		self.timepoints = []
		self.concs = []
		self.num_molecules = []
		self.ids_label = []


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
			- number (numpy[int]): Time series of the numbers of molecules in labels (3D array; [Timepoints, species, ids_spine]).
			- conc in uM (numpy[float]): Time series of the concentrations of molecules in labels (unit: uM) (3D array; [Timepoints, species, ids_spine]).
		"""
		if not isinstance(filename, str):
			raise ValueError('filename is not str.')
		elif (self.timepoints == []) or (self.concs == []) or (self.numbers == []) or (self.ids_label == []):
			raise ValueError('No labeled concs. exec_label_volume or exec_label_file.')

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
			f['ids_label']  = self.ids_label
			f['species']    = species
			f['numbers']    = data_numbers
			f['concs in uM'] = data_concs


	def exec_volume(self, label_volume):
		"""Get time series of moleuclar numbers/concentrations within labeled volumes from LM simulation result.

		Args:
			label_volume (numpy[int/bool]): Label volume (3D array). id=0 will be ignored. label_filename or label_volume must be specified.

		Returns:
			(obj): Instance variables containing:

			- timepoints (numpy[float]): Timepoints (s)
			- concs (dict): Time series of molecular concentrations (number per labeled volume, in the unit of uM). The dict container has {'species1': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], 'Species2': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], ... }.
			- numbers (dict): Time series of the numbers of molecules of the specified molecular species. The dict container has {'species1': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], 'Species2': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], ... }.
			- ids_spine (numpy[int]): Labels in the volume. The numpy array has [label1, label2, ...].
		"""
		self._exec(label_volume)


	def exec_file(self, label_filename):

		"""Get time series of moleuclar numbers/concentrations within labeled volumes from LM simulation result.

		Args:
			label_filename (str): Filename of label. label_filename or label_volume must be specified.

		Returns:
			(obj): Instance variables containing:

			- timepoints (numpy[float]): Timepoints (s)
			- concs (dict): Time series of molecular concentrations (number per labeled volume, in the unit of uM). The dict container has {'species1': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], 'Species2': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], ... }.
			- numbers (dict): Time series of the numbers of molecules of the specified molecular species. The dict container has {'species1': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], 'Species2': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], ... }.
			- ids_spine (numpy[int]): Labels in the volume. The numpy array has [label1, label2, ...].
		"""
		if isinstance(label_filename, str) :
			with h5py.File(label_filename,'r') as f:
				label_volume = f['label volume'][()]
		else:
			print('label_filename must be str.')
			return False

		self._exec(label_volume)

	def _exec(self, label_volume):

		# Check arguments
		if not isinstance(label_volume, np.ndarray) or (label_volume.ndim != 3) or (label_volume.dtype not in NUMPY_INTEGERS):
			raise ValueError('label_volume must be a integer 3D np.ndarray.')

		ids_label, nums_label_voxels = np.unique(label_volume, return_counts=True)
		ids_label          = ids_label[1:] #  0 was removed.
		nums_label_voxels  = nums_label_voxels[1:]

		targ_spine_labels = []
		labels_flat = np.ravel(label_volume)
		for id_label in ids_label:
		    targ_spine_label = np.where( (labels_flat == id_label) )
		    targ_spine_labels.append( targ_spine_label )

		# Time frames
		with h5py.File(self.lm_filename,'r') as file:
			timepoints = file['Simulations']['0000001']['LatticeTimes'][()]
			frames = [key for key in file['Simulations']['0000001']['Lattice'].keys()]
		frames.sort()

		# Monitor messages
		if self.monitor_species != None:
		    print('file :', lm_filename)
		    print('Label ids : ', ids_label)
		    print('Species    : ', ', '.join(list(S.keys())))


		# Get molecules in labeled volumes at each timepoint
		numbers = {}
		for Targ in self.s.keys():
		    numbers[Targ] = []

		for f in frames:
			#
			num_molecules_time_i = {}
			for Targ in self.s.keys():
			    num_molecules_time_i[Targ] = [0 for i in range(max(ids_label)+1)]

			#
			with h5py.File(self.lm_filename,'r') as file:
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
		self.ids_label = ids_label


