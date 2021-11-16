import sys, os
import h5py

class LoadLabelVolume():

	"""A utility class to load label info from a label volumes file.

	Args:
	    h5_filename (str): HDF5 filename

	Returns:
		(pyLD.LoadLabeledVolumes): LoadLabeledVolumes object that contains the following instance variables:

		- 'label_volume' (numpy[int]): Loaded label volume (3D array) 
		- 'label_ids' (numpy[int]): Ids of new labels (1D array)
		- 'ref_volume' (numpy[bool]): Reference domain volume  (3D array)
	"""
	def __init__(self, h5_filename):
		if not isinstance(h5_filename, str):
			raise ValueError('h5_filename must be str.')
		self._load(h5_filename)


	def _load(self, h5_filename):
		with h5py.File(h5_filename, 'r') as f:
			self.label_volume = f['label volume'][()]
			self.label_ids    = f['label ids'][()]
			self.ref_volume   = f['ref volume'][()]
