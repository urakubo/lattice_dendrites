from __future__ import print_function
from __future__ import division
import numpy as np
import os
import random


NUMPY_INTEGERS = [ np.int8, np.int16, np.int32, np.int64, \
	np.uint8, np.uint16, np.uint32, np.uint64]



class BuildAnyShape:
	"""Build a buildAnyShape object.

	sim.siteTypes['domain name'] may be different from volume_id.

	Args:
		sim (obj): RDMESimulation object
		volume (numpy[int]): 3D array that specifies volume_ids
		domains (dict): {'domain name', volume_id}
		surfaces (dict): {'surface_name', numpy[float]} : Surface areas in voxel space (3D array)

	Returns:
		(pyLD.buildAnyShape): buildAnyShape object
	"""

	def __init__(self, sim, volume, domains, surfaces = {}):

		self._check_arguments(sim, volume, domains, surfaces)
		self.sim = sim

		# Add regions and rename regions of the target volume
		volume_mod = volume.astype(np.int) * 0
		for domain_name in domains:
		    sim.addRegion( domain_name )
		    volume_mod += (volume == domains[domain_name]).astype(np.int) * sim.siteTypes[domain_name]

		# Check volume size
		x_lsize = sim.lattice.getXSize()
		y_lsize = sim.lattice.getYSize()
		z_lsize = sim.lattice.getZSize()
		print('Lm lattice size   (x, y, z): ', x_lsize, y_lsize, z_lsize )
		x_vsize  = volume.shape[0]
		y_vsize  = volume.shape[1]
		z_vsize  = volume.shape[2]
		print('Input volume size (x, y, z): ', x_vsize, y_vsize, z_vsize )
		xnum   = min([x_lsize, x_vsize])
		ynum   = min([y_lsize, y_vsize])
		znum   = min([z_lsize, z_vsize])

		# Set volume
		# Can be accelerated using serialization.
		for x in range(xnum):
		    for y in range(ynum):
		        for z in range(znum):
		            sim.lattice.setSiteType(x, y, z, int(volume_mod[x,y,z]))
		sim.hasBeenDiscretized = True

		# Register domain locations.
		## dict.fromkeys(domains, []) refers an identical empty list []. 
		self.locs   = dict.fromkeys(domains)
		for domain_name in domains:
		    self.locs[domain_name] = []

		for x in range(xnum):
		    for y in range(ynum):
		        for z in range(znum):
		            for domain_name in domains:
		                if (volume_mod[x,y,z] == sim.siteTypes[domain_name]):
		                    self.locs[domain_name].append((x,y,z))


		self.num_voxels = dict.fromkeys(domains)
		for domain_name in domains:
		    self.num_voxels[domain_name] = len(self.locs[domain_name])


		# Extract the regions of membranes and PSDs
		
		self.surf_voxel_locs = {}
		self.surf_voxel_prob = {}
		
		for k, v in surfaces.items():
			self.surf_voxel_locs[k] = np.nonzero(v > 0)
			self.surf_voxel_prob[k] = v[self.surf_voxel_locs[k][0], self.surf_voxel_locs[k][1], self.surf_voxel_locs[k][2] ]

		# print('self.memb_voxel_locs[0].shape[0]: ', self.memb_voxel_locs[0].shape[0])
		#print("len(locs['default'])   : ", len(self.locs['default']))
		#print("len(locs['cytoplasm']) : ", len(self.locs['cytoplasm']))
		#print("len(locs['psd'])       : ", len(self.locs['psd']))                           

	def add_solute_molecules(self, molecular_name, molecular_number, domain_name):
		"""Add solute molecules.

		Args:
			molecular_name (str): Molecular name
			molecular_number (int): Number of molecules
			domain_name (str): Domain name

		Returns:
			(bool): The return value. True for success, False otherwise.
		"""

		particleNum=self.sim.particleMap[molecular_name]
		for i in random.sample(self.locs[domain_name], molecular_number):
		    self.sim.lattice.addParticle(i[0], i[1], i[2], int(particleNum))
		self.sim.customAddedParticleList.append((molecular_name, molecular_number))
		return True


	def add_surface_molecules(self, molecular_name, density, surface_name):
		"""Add membrane molecules.

		Args:
			molecular_name (str): Molecular name
			density (float): Number density (/um2?)
			surface_name (str): Surface name

		Returns:
			(bool): The return value. True for success, False otherwise.
		"""

		particleNum=self.sim.particleMap[molecular_name]
		molecular_numbers = np.random.binomial(density, self.surf_voxel_prob[surface_name])
		# print('np.sum(molecular_numbers): ', np.sum(molecular_numbers))
		locs = self.surf_voxel_locs[surface_name]
		for x, y, z, num in zip(locs[0], locs[1], locs[2], molecular_numbers):
		    # print('x,y,z, num: ', x,y,z, num)
		    for i in range(num):
		        self.sim.lattice.addParticle(int(x), int(y), int(z), int(particleNum))
		return True


	def _check_arguments(self, sim, volume, domains, surfaces):


		try:
		    import pySTDLM.StandardReactionSystems
		    if type(sim) != pySTDLM.StandardReactionSystems.RDMESimulation:
		        raise ValueError('volume must be a integer 3D np.ndarray.')

		except ImportError:
		    pass

		if not isinstance(volume, np.ndarray) or (volume.ndim != 3) or (volume.dtype not in NUMPY_INTEGERS):
			raise ValueError('volume must be a integer 3D np.ndarray.')
		elif not isinstance(domains, dict) :
			raise ValueError('domains must be a dict.')
		elif not isinstance(surfaces, dict) :
			raise ValueError('surfaces must be a dict.')
		else:
			return True

