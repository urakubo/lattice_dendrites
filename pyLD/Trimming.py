class Trimming():
	"""Trim unused areas from a volume to obtain a minimal bounding box.

	Args:
	    reference_volume (numpy[int/bool]): Reference volume to calcutate a minimal bounding box
	"""
	def __init__(self, reference_volume):
		self.rmin, self.rmax, self.cmin, self.cmax, self.zmin, self.zmax = self._bbox_3D(reference_volume)

	def trim(self, volume):
		"""Execute the trimming of a specified volume.

		Args:
	    	volume (numpy[int/bool]): Target volume

		Returns:
			(numpy[int/bool]): Rotated volume
		"""
		volume = volume[self.rmin:self.rmax, self.cmin:self.cmax, self.zmin:self.zmax]
		return target_volume

	def _bbox_3D(self, vol):
		r = np.any(vol, axis=(1, 2))
		c = np.any(vol, axis=(0, 2))
		z = np.any(vol, axis=(0, 1))
		rmin, rmax = np.where(r)[0][[0, -1]]
		cmin, cmax = np.where(c)[0][[0, -1]]
		zmin, zmax = np.where(z)[0][[0, -1]]

		return rmin, rmax, cmin, cmax, zmin, zmax
