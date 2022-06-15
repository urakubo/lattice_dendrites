import numpy as np
import h5py
import os, sys
from .utils import get_species_names, get_spacing, num_to_uM


class ConnectTotalConcs:
	"""Obtain time developments of the total numbers/concentrations of target molecules.
	
	Args:
		label_conc_filenames (str / list[str] / tuple[str]): Simulated lm files.
		domain_names (str / list[str] / tuple[str]): Target domain names. They are summed if multiple domains are specified.

	Returns:
		(pyLD.ConnectTotalConcs): ConnectTotalConcs object that has the follwing instances:

		- timepoints (numpy[float]): Timepoints (s)
		- concs (numpy[float]): Time series of the concentrations of molecules in labels (unit: uM) (3D array; [Timepoints, species, ids_label]).
		- numbers (numpy[float]): Time series of the numbers of molecules in labels (3D array; [Timepoints, species, ids_label]).
		- ids_label (numpy[int]): Labels. It has [label1, label2, ...].
	"""

	def  __init__(self, lm_files, domain_names):
		# Check arguments
		if isinstance(lm_files, str):
		    lm_files = [lm_files]
		elif isinstance(lm_files, list) | isinstance(lm_files, tuple) :
		    pass
		else:
		    raise ValueError('lm_files must be str, list, or tuple.')

		if isinstance(domain_names, str):
		    domain_names = [domain_names]
		elif isinstance(domain_names, list) | isinstance(domain_names, tuple) :
		    pass
		else:
		    raise ValueError('domain_names must be str, list[str], or tuple[str].')


		# Timepoints
		for i, fname in enumerate(lm_files):
		    with h5py.File(fname,'r') as f:
		        t = f['Simulations']['0000001']['LatticeTimes'][()]
		        n = f['Simulations']['0000001']['SpeciesCounts'][()]
		    if i == 0:
		        self.timepoints = t
		        self.numbers    = n
		    else:
		        self.timepoints = np.hstack( (self.timepoints, t[1:]+self.timepoints[-1]) )
		        self.numbers    = np.concatenate([self.numbers , n[1:,:] ], axis=0)

		# Concs

		with h5py.File(lm_files[0],'r') as f:
		    spacing  = f['Model']['Diffusion'].attrs['latticeSpacing']
		    volume   = f['Model']['Diffusion']['LatticeSites'][()]
		    siteType = f['Parameters'].attrs['siteTypeNames'].decode().split(',')
		    self.species = f['Parameters'].attrs['speciesNames'].decode().split(',')
		num_voxels  = 0
		for name in domain_names:
			num_voxels += np.count_nonzero(volume == siteType.index(name))
		self.concs = num_to_uM(self.numbers, num_voxels, spacing)


	def get_concs(self, species):
		"""Get time developments of the concentration(s) of specified specie(s).

		Args:
			species (str / list[str] / tuple[str]): Target molecular species. They are summed if multiple species are specified.

		Returns:
			(numpy[float]): Time developments of concentrations (1D array)
		"""

		species_ids = self._check_arguments(species)
		return np.sum(self.concs[:, species_ids], axis = 1)


	def get_numbers(self, species):
		"""Get time developments of the number(s) of specified specie(s).

		Args:
			species (str / list[str] / tuple[str]): Target molecular species. They are summed if multiple species are specified.

		Returns:
			(numpy[int]): Time developments of number (1D array)
		"""

		species_ids = self._check_arguments(species)
		return np.sum(self.numbers[:, species_ids], axis = 1)


	def _check_arguments(self, species):
		# Check arguments
		if isinstance(species, str):
		    species = [species]
		elif isinstance(species, list) | isinstance(species, tuple) :
			pass
		else:
			raise ValueError('species must be str, list[str], or tuple[str].')

		ids = [i for i, key in enumerate(self.species) if key in species ]
		return ids


