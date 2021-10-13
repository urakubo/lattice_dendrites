from __future__ import print_function
from __future__ import division

import numpy as np
import cv2
# import matplotlib.pyplot as plt

class RotateVolume():
	"""Rotate a volume to a minimal bounding box using cv2.minAreaRect.
	Domains are transformed in a minimal bounding box of a reference volume in the X-Y space.

	Args:
	    reference_volume (numpy[int/bool]): Reference volume to calcutate a minimal bounding box
	    fixed_axis (int): Fixed axis in rotation (x:0, y:1, z:2)
	"""
	def __init__(self, reference_volume, fixed_axis = 0):

		self.d  = 2
		self._calc(reference_volume, fixed_axis)

	def _calc(self, volume, fixed_axis):
		"""Calculate a rotation matrix.
		_calc is automatically called from the initial definition.
		Users can redefine the rotation matrix.

		Args:
	    	volume (numpy[int/bool]): Reference volume to calcutate a minimal bounding box
		    fixed_axis (int): Fixed axis in rotation (x:0, y:1, z:2)
		"""

		self.fixed_axis = fixed_axis
		volume = volume.swapaxes(0, self.fixed_axis)

		summed_image = np.sum((volume > 0).astype(np.int), axis=0)
		summed_image = (summed_image > 0).astype(np.uint8)
		if np.sum(summed_image) == 0:
			print('No domain.')
			return False

		self.cols, self.rows   = summed_image.shape[0], summed_image.shape[1]
		_contours,_ = cv2.findContours(summed_image, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
		contours = _contours[0]
		for i in range(1,len(_contours)):
			contours = np.concatenate([contours, _contours[i]], 0)
		contours = contours[:,0,:]

#		plt.scatter(contours[:,0], contours[:,1])
#		plt.show()
		rect = cv2.minAreaRect(contours)
		#print('rect: ', rect)
		#(center(x, y), (width, height), angle of rotation) 
		
		# rotation matrix
		theta   = rect[2] * 3.14159/180# [0, -90 degree)
		self.M = np.array([[np.cos(-theta), -np.sin(-theta)],[np.sin(-theta), np.cos(-theta)]])
		# Bounding box rotation
		box_before = cv2.boxPoints(rect)
		box_before = np.array(box_before).T
		box_after  = self.M @ box_before
		self.xmin_box_after = np.min(box_after[0,:]).astype(int)
		self.xmax_box_after = np.max(box_after[0,:]).astype(int)
		self.ymin_box_after = np.min(box_after[1,:]).astype(int)
		self.ymax_box_after = np.max(box_after[1,:]).astype(int)
		self.xwidth_box_after = self.xmax_box_after - self.xmin_box_after
		self.ywidth_box_after = self.ymax_box_after - self.ymin_box_after

		self.M = np.array([[np.cos(-theta), -np.sin(-theta), -self.xmin_box_after], \
				[np.sin(-theta), np.cos(-theta), -self.ymin_box_after]])

#		self.M = cv2.getRotationMatrix2D( (self.xmin_box_after, self.ymin_box_after), theta, 1)


	def rotate(self, volume):
		"""Rotate a specified volume.

		Args:
	    	volume (numpy[int/bool]): Target volume

		Returns:
			(numpy[int/bool]): Rotated volume
		"""
		volume_dtype  = volume.dtype
		volume_transp = volume.swapaxes(0, self.fixed_axis)

		volume_transp_rotated = np.zeros((volume_transp.shape[0]+2*self.d, \
					self.ywidth_box_after+2*self.d, \
					self.xwidth_box_after+2*self.d), dtype=volume_dtype)

		for i in range(volume_transp.shape[0]):
			slice     = volume_transp[i,:,:].astype(np.int)
			slice_rot = cv2.warpAffine( slice, \
					self.M, (self.cols+self.xwidth_box_after, \
					self.rows + self.ywidth_box_after), \
					flags=cv2.INTER_NEAREST )
			volume_transp_rotated[i+self.d, self.d:-self.d, self.d:-self.d] = \
					slice_rot[0: self.ywidth_box_after, 0: self.xwidth_box_after ].astype(volume_dtype)


		return volume_transp_rotated

