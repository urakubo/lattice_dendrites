#
# AquireVolumeFromDxfFile
#


import sys, os, errno
import glob
import numpy as np
import cv2
import pathlib
from natsort import natsorted
import copy
import matplotlib.pyplot as plt
 
class createVolumeFromReconstructFile():

	def __init__(self, path_to_reconstruct_dxf_files, xypitch, zpitch):

		if zpitch % xypitch != 0:
			print('zpitch must be a multiple of xypitch!')
			return False

		self.rotation_matrix = []
		self.xypitch = xypitch
		self.zmult = int(zpitch / xypitch)

		self.d  = 2
		self.dz = 1

		p = pathlib.Path(path_to_reconstruct_dxf_files)
		if p.is_absolute() == False:
			self.path_to_dxf = p.resolve()
		else:
			self.path_to_dxf = path_to_reconstruct_dxf_files

		self.dxf_files = natsorted( glob.glob( os.path.join(self.path_to_dxf, "*.dxf") ) )
		if self.dxf_files == []:
			print('No reconstruct dxf files!')
			return False

		# print(self.dxf_files)


	def calc_bounding_box(self, target_domain):

		if isinstance(target_domain, list) or isinstance(target_domain, tuple):
			domains = target_domain
		elif isinstance(target_domain, str):
			domains = [target_domain]
		else:
			print('Target domains must be str, tuple, or list!, but:')
			print(target_domain)
			return False

		ids = []
		vertices = []
		for i, dxf_file in enumerate(self.dxf_files):
			for domain in domains:
				vvs = self._load_vertices_xy(dxf_file, domain)
				if vvs != []:
					ids.append(i)
					for vv in vvs:
						vertices.extend(vv)
		self.max_id_z = max(ids)
		self.min_id_z = min(ids)

		if vertices == []:
			print('No target domain: ', target_domain)
			return False

		vertices = (np.array(vertices)/self.xypitch).astype('int')
		self.rect = cv2.minAreaRect(vertices)
		box = cv2.boxPoints(self.rect)

#       (center(x, y), (width, height), angle of rotation) = cv2.minAreaRect(points)
#		print(self.rect)

		# rotation matrix
		angle = self.rect[2]
		transformation_matrix = cv2.getRotationMatrix2D((0,0) ,angle, 1)
		self.rotation_matrix  = np.linalg.inv(transformation_matrix[:2,:2])

		# Transformation of boundingbox
		translated = box @ self.rotation_matrix
		self.minmax_col = [np.min(translated[:,0]), np.max(translated[:,0])]
		self.minmax_row = [np.min(translated[:,1]), np.max(translated[:,1])]

		# Set volume size
		self.col_volume = (self.minmax_col[1] - self.minmax_col[0] + self.d * 2).astype(int)
		self.row_volume = (self.minmax_row[1] - self.minmax_row[0] + self.d * 2).astype(int)
		self.z_volume   = (self.max_id_z - self.min_id_z + 1)*self.zmult + self.dz * 2


	def create(self, target_domain):

		if self.rotation_matrix == []:
			print('Calculate bounding_box beforehand.')
			return False

		if isinstance(target_domain, list) or isinstance(target_domain, tuple):
			domains = target_domain
		elif isinstance(target_domain, str):
			domains = [target_domain]
		else:
			print('Target domains must be str, tuple, or list!, but:')
			print(target_domain)
			return False


		### Why ...???
		volume    = np.zeros((self.row_volume, self.col_volume, self.z_volume), dtype='bool')
		slice_ref = np.zeros((self.row_volume, self.col_volume), dtype='uint8')

		for i in range(self.min_id_z, self.max_id_z):

			slice =	copy.deepcopy(slice_ref)

			pp = []
			for domain in domains:
				vvs = self._load_vertices_xy(self.dxf_files[i], domain)
				if vvs == []:
					continue
				# print(dxf_file)
				for vs in vvs:
					p = (np.array(vs)/self.xypitch)
					p = p @ self.rotation_matrix -[self.minmax_col[0], self.minmax_row[0]]
					p = p.astype(np.int32) + [self.d, self.d]
					pp.append(p)

			slice = cv2.fillPoly(slice, pp, 1)
			for j in range(self.zmult):
				volume[:,:,(i-self.min_id_z)*self.zmult + j + self.dz] = slice.astype(np.bool)

		return volume


#
#
	def _load_vertices_xy(self, dxf_file, target_domain):
		
		with open(dxf_file, 'r') as f:
			textlist = f.readlines()

		# IDs & NUM of POLYLINE
		idpolys      = [i for i, text in enumerate(textlist) if 'POLYLINE'in text]
		namepolys    = [textlist[i+8] for i in idpolys]
		# End of POLYLINE
		idpolyends  = [i for i, text in enumerate(textlist) if 'SEQEND'in  text]
		# IDs of Vertices
		idverts     = [i for i, text in enumerate(textlist) if 'VERTEX' in text]

		polylines = []
		for idpoly, namepoly, idpolyend in zip(idpolys, namepolys, idpolyends):
			if namepoly != target_domain+'\n':
				continue
			else:
				ids = [id for id in idverts if (id > idpoly) and (id < idpolyend)]
				polyline = [[float(textlist[i+2]), float(textlist[i+4])] for i in ids]
				polylines.append(polyline)

		return polylines



			# tmp = cv2.fillConvexPoly(slice_ref, p, 25)


			#p = p.reshape((-1,1,2)).astype(np.int32)
			# campus = cv2.drawContours(campus,pp,0,255,3)
			# print(pp)


class rotateVolume():

	def __init__(self, volume, fixed_axis = 0):

		self.d  = 2
		self.calc(volume, fixed_axis)

	def calc(self, volume, fixed_axis):

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


	def exec(self, volume):

		volume_dtype  = volume.dtype
		volume_transp = volume.swapaxes(0, self.fixed_axis)


		volume_transp_rotated = np.zeros((volume_transp.shape[0]+2*self.d, \
					self.ywidth_box_after+2*self.d, \
					self.xwidth_box_after+2*self.d), dtype=volume_dtype)

		for i in range(volume_transp.shape[0]):
			slice     = volume_transp[i,:,:].astype(np.int)
			slice_rot = cv2.warpAffine( slice, self.M, (self.cols+self.xwidth_box_after, self.rows + self.ywidth_box_after) , flags=cv2.INTER_NEAREST )
			volume_transp_rotated[i+self.d, self.d:-self.d, self.d:-self.d] = slice_rot[0: self.ywidth_box_after, 0: self.xwidth_box_after ].astype(volume_dtype)

			"""
			if i == 100:
				fig, (ax1, ax2, ax3) = plt.subplots(1,3)
				ax1.imshow(slice, cmap="gray")
				ax2.imshow(slice_rot, cmap="gray")
				ax2.plot([0, self.xwidth_box_after, self.xwidth_box_after, 0,  0], \
						[0, 0, self.ywidth_box_after, self.ywidth_box_after, 0],  'r-')

				ax3.imshow(volume_transp_rotated[i,:,:], cmap="gray")
				plt.show()
			"""

		return volume_transp_rotated


#
#
#
import collections

def flatten(l):
	for el in l:
		if isinstance(el, collections.abc.Iterable) and not isinstance(el, (str, bytes)):
			yield from flatten(el)
		else:
			yield el


