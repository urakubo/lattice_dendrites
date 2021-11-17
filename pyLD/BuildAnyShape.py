from __future__ import print_function
from __future__ import division
import numpy as np
import os
import random
from .utils import uM_to_num
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
		(pyLD.BuildAnyShape): BuildAnyShape object
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
                
		# unit converter
		NA           = 6.022e23
		self.spacing = sim.latticeSpacing
		experiment_volume_in_L = self.spacing * self.spacing * self.spacing * \
                                         1000 * x_lsize * y_lsize * z_lsize
		self.per_uM = 1/(NA*1e-6*experiment_volume_in_L)


                
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


                
	def number_per_1uM(self, domain_name ):
		num_voxels = self.num_voxels[domain_name]
		spacing_in_m = self.spacing
		conc_in_uM = 1
		number_per_1uM = uM_to_num(conc_in_uM, num_voxels, spacing_in_m)
		return number_per_1uM

                
	def define_species(self, molecular_names):
		# Check arguments
		if isinstance(molecular_names, str):
		    molecular_names = [molecular_names]
		elif not isinstance(molecular_names, list) and not isinstance(molecular_names, tuple) :
		    raise ValueError('lm_files must be str, list, or tuple.')
		else:
		    self.sim.defineSpecies(molecular_names)

	def modify_region(self, domain_name ):
		#
		return self.sim.modifyRegion( domain_name )

	def modifyRegion(self, domain_name ):
		#
		return self.sim.modifyRegion( domain_name )

	def add_reaction(self, reactant, product, rate, domain_name ):
		#
		self.sim.modifyRegion( domain_name ).addReaction(reactant=reactant, product=product, rate=rate)

	def reac_oneway_uM(self, reactant, product, rate, domain_name ):
		if len(reactant) == 2:
			self.sim.modifyRegion(domain_name).addReaction(reactant=reactant, product=product, rate=rate*self.per_uM)
		else:
			self.sim.modifyRegion(domain_name).addReaction(reactant=reactant, product=product, rate=rate)

	def reac21_uM(self, reac1, reac2, prod, kf, kb, domain_name ):
		self.sim.modifyRegion(domain_name)\
				.addReaction(reactant=(reac1, reac2), product=prod, rate=kf*self.per_uM)
				.addReaction(reactant=prod, product=(reac1, reac2), rate=kb)

	def reac12_uM(self, reac, prod1, prod2, kon, koff, domain_name ):
		self.sim.modifyRegion(domain_name)\
				.addReaction(reactant=reac, product=(prod1, prod2), rate=kf)
				.addReaction(reactant=(prod1, prod2), product=reac, rate=kb*self.per_uM)

	def set_diffusion_rate(self, molecular_name, rate, domain_name ):
		#
		self.sim.modifyRegion( domain_name ).setDiffusionRate(molecular_name, rate=rate)

	def add_solute_molecules_uM(self, molecular_name, number, domain_name):
		#
		number_per_1uM = self.number_per_1uM(domain_name)
		self.add_solute_molecules( molecular_name, int(number*number_per_1uM), domain_name )

	def add_solute_molecules(self, molecular_name, number, domain_name):
		"""Add solute molecules.

		Args:
			molecular_name (str): Molecular name
			number (int): Number of molecules
			domain_name (str): Domain name

		Returns:
			(bool): The return value. True for success, False otherwise.
		"""
		particleNum=self.sim.particleMap[molecular_name]
		for i in random.sample(self.locs[domain_name], number):
		    self.sim.lattice.addParticle(i[0], i[1], i[2], int(particleNum))
		self.sim.customAddedParticleList.append((molecular_name, number))
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
			raise ValueError('domains must be dict.')
		elif not isinstance(surfaces, dict) :
			raise ValueError('surfaces must be dict.')
		else:
			return True

