#
# AquireVolumeFromDxfFile
#
from __future__ import print_function
from __future__ import division

import sys, os, errno
import glob
import numpy as np
import cv2
import pathlib
import copy
import re

#import matplotlib.pyplot as plt
# from natsort import natsorted

from skimage.transform import rescale

 
class createVolumeFromReconstruct():
	"""Create a volume from a set of Reconstruct dxf files.
	Multiple domains can be converted into a volume.
	Domains are transformed into a minimal bounding box of the X-Y space.

	Args:
		path_to_reconstruct_dxf_files (str): Relative or absolute path to dxf files
		xypitch (float): XYZ-pitch per voxel
		zpitch (float): Z-slice thickness
		reference_domains (str/list[str]/tuple[str]): Reference domain names to calcutate a minimal bounding box

	Returns:
		(pyLD.createVolumeFromReconstruct): createVolumeFromReconstruct object
	"""

	def __init__(self, path_to_reconstruct_dxf_files, xypitch, zpitch, reference_domains):

		self.rotation_matrix = []
		self.xypitch = xypitch
		self.zmult = zpitch / xypitch

		self.d  = 2
		self.dz = 1

		p = pathlib.Path(path_to_reconstruct_dxf_files)
		if p.is_absolute() == False:
			self.path_to_dxf = p.resolve()
		else:
			self.path_to_dxf = path_to_reconstruct_dxf_files

#		self.dxf_files = natsorted( glob.glob( os.path.join(self.path_to_dxf, "*.dxf") ) )
		self.dxf_files = glob.glob( os.path.join(self.path_to_dxf, "*.dxf") )
		self.dxf_files = sorted(self.dxf_files, key=lambda s: int(re.findall(r'\d+', s)[-1]))
		if self.dxf_files == []:
			print('No reconstruct dxf files!')
			return False

		# print(self.dxf_files)
		self._calc_bounding_box(reference_domains)


	def _calc_bounding_box(self, reference_domains):
		"""Calculate a transformation matrix for a domain to fit a minimal bounding box.
		_calc_bounding_box is automatically called from the initialization.
		Users can redefine the rotation matrix.

		Args:
			reference_domains (str/list[str]/tuple[str]): Reference domain names to calcutate a minimal bounding box
		Returns:
			(pyLD.createVolumeFromReconstruct): createVolumeFromReconstruct object
		"""

		if isinstance(reference_domains, list) or isinstance(reference_domains, tuple):
			domains = reference_domains
		elif isinstance(reference_domains, str):
			domains = [reference_domains]
		else:
			print('Target domains must be str, tuple, or list!, but:')
			print(reference_domains)
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
			print('No reference domain(s): ', reference_domains)
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
		self.z_volume   = self.max_id_z - self.min_id_z + self.dz * 2


	def create(self, target_domain):
		"""Create a volume of a target domain.

		Args:
			target_domain (str): Target domain name
		Returns:
			(numpy[bool]): Volume
		"""

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
		volume    = np.zeros((self.row_volume, self.col_volume, self.z_volume), dtype='uint8')
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

			slice = cv2.fillPoly(slice, pp, color = 255)
			volume[:,:, i-self.min_id_z + self.dz] = slice.astype( np.uint8 )

		#
		# Expansion in the Z axis.
		#
		volume = rescale(volume, scale=(1, 1, self.zmult), order=3, preserve_range=True, multichannel=False) # order 1 = biliner,  order 3 = bicubic
		volume = (volume > 127).astype(np.bool)

		return volume


	def _load_vertices_xy(self, dxf_file, target_domain):
		"""Load vertices of target_domain from a dxf_file.

		Args:
			dxf_file (str): Dxf file name
			target_domain (str): Target domain name
		Returns:
			(list[list[float,float]]): List of vertices
		"""

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


