import numpy as np
import h5py
import os, sys
from .utils import get_species_names, get_spacing, num_to_uM

NUMPY_INTEGERS = [ np.int8, np.int16, np.int32, np.int64, \
	np.uint8, np.uint16, np.uint32, np.uint64]


class ConnectLabeledConcs:
	"""Connect the time developments of molecular numbers/concentrations in labeled conc files.
	
	Args:
		label_conc_filenames (str / list[str] / tuple[str]): Generated labeled conc files from "get_labeled_concs".

	Returns:
		(pyLD.ConnectLabeledConcs): ConnectLabeledConcs object that includes the follwing instances:

		- timepoints (numpy[float]): Timepoints (s)
		- concs (numpy[float]): Time series of the concentrations of molecules in labels (unit: uM) (3D array; [Timepoints, species, ids_label]).
		- numbers (numpy[float]): Time series of the numbers of molecules in labels (3D array; [Timepoints, species, ids_label]).
		- ids_label (numpy[int]): Labels. It has [label1, label2, ...].
	"""

	def  __init__(self, label_conc_filenames):
		# Check arguments
		if isinstance(label_conc_filenames, str):
		    label_conc_filenames = [label_conc_filenames]
		elif isinstance(label_conc_filenames, list) | isinstance(label_conc_filenames, tuple) :
			pass
		else:
			raise ValueError('label_conc_filenames must be str, list, or tuple.')

		# Species names, ids_spine
		self.filenames = label_conc_filenames
		with h5py.File(self.filenames[0],'r') as f:
			self.species   = list(f['species'][()])
			self.ids_label = f['ids_label'][()]

		# Timepoints
		for i, fname in enumerate(self.filenames):
		    with h5py.File(fname,'r') as f:
		        t = f['timepoints'][()]
		        c = f['concs in uM'][()]
		        n = f['numbers'][()]
		    if i == 0:
		        self.timepoints = t
		        self.concs      = c
		        self.numbers    = n
		    else:
		        self.timepoints = np.hstack( (self.timepoints, t[1:]+self.timepoints[-1]) )
		        self.concs      = np.concatenate([self.concs   , c[1:,:,:] ], axis=0)
		        self.numbers    = np.concatenate([self.numbers , n[1:,:,:] ], axis=0)


	def get_concs(self, species=None, ids_label=None):
		"""Get time developments of the concentration(s) of specified specie(s) within label(s).

		Args:
			species (None / str / list[str] / tuple[str]): Molecular species. They are summed if multiple species are specified, and unsummed if None is specified.
			ids_label (None / int / list[int] / tuple[int]): Label ids. They are summed if multiple labels are specified, and unsummed if None is specified.

		Returns:
			(numpy[float]): Time developments of concentrations (1D/2D/3D array)
		"""

		species_id, label_id = self._check_arguments(species, ids_label)

		if (species != None) and (ids_label != None):
			concs = self.concs[:,species_id,label_id]
			if concs.ndim == 3:
				concs = np.sum(concs, axis=(1,2))
			elif concs.ndim == 2:
				concs = np.sum(concs, axis=1)
		elif (species != None):
			concs = np.sum(self.concs[:,species_id,:], axis=(1))
		elif (ids_label != None):
			concs = np.sum(self.concs[:,:,label_id], axis=(2))
		else:
			concs = self.concs

		return concs


	def get_numbers(self, species=None, ids_label=None):
		"""Get time developments of the number(s) of specified specie(s) within label(s).

		Args:
			species (None / str / list[str] / tuple[str]): Molecular species. They are summed if multiple species are specified, and unsummed if None is specified.
			ids_label (None / int / list[int] / tuple[int]): Label ids. They are summed if multiple labels are specified, and unsummed if None is specified.

		Returns:
			(numpy[int]): Time developments of number (1D/2D/3D array)
		"""

		species_id, label_id = self._check_arguments(species, ids_label)

		if (species != None) and (ids_label != None):
			numbers = self.numbers[:,species_id,label_id]
			if numbers.ndim == 3:
				numbers = np.sum(numbers, axis=(1,2))
			elif numbers.ndim == 2:
				numbers = np.sum(numbers, axis=1)
		elif (species != None):
			numbers = np.sum(self.numbers[:,species_id,:], axis=(1))
		elif (ids_label != None):
			numbers = np.sum(self.numbers[:,:,label_id], axis=(2))
		else:
			numbers = self.numbers

		return numbers


	def _check_arguments(self, species, ids_label):

		# Check arguments
		if isinstance(species, str):
		    species = [species]
		elif (species == None) or isinstance(species, list) or isinstance(species, tuple) :
			pass
		else:
			raise ValueError('species must be None, str, list, or tuple.')

		if isinstance(ids_label, int):
		    ids_label = [ids_label]
		elif (ids_label == None) or isinstance(ids_label, list) or isinstance(ids_label, tuple) :
			pass
		else:
			raise ValueError('ids_label must be None, int, list, or tuple.')


		if (species != None):
			species_id = [id for id, s in enumerate(self.species) if s in species]
		else:
			species_id = None

		if (ids_label != None):
			label_id = [id for id, i in enumerate(self.ids_label) if i in ids_label ]
		else:
			label_id = None

		# print('species_id: ', species_id)
		# print('label_id  : ', label_id)
		return species_id, label_id


