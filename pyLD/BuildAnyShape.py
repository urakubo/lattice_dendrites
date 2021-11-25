from __future__ import print_function
from __future__ import division
import numpy as np
import os
import random
from .utils import uM_to_num
NUMPY_INTEGERS = [ np.int8, np.int16, np.int32, np.int64, \
	np.uint8, np.uint16, np.uint32, np.uint64]

try:
	import pyLM.RDME
	import pyLM.units
except ImportError:
    pass


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

	def __init__(self, volume, domains, spacing_in_m, surfaces = {}):

		self._check_arguments(volume, domains, spacing_in_m, surfaces)

		nx, ny, nz = volume.shape
		self.sim = pyLM.RDME.RDMESimulation(
			dimensions= (nx * spacing_in_m, ny * spacing_in_m, nz * spacing_in_m),\
			spacing=spacing_in_m)

		self.spacing_in_m = spacing_in_m
		self.nx = nx
		self.ny = ny
		self.nz = nz

		# Add regions and rename regions of the target volume
		volume_mod = volume.astype(np.int) * 0
		for domain_name in domains:
		    self.sim.addRegion( domain_name )
		    volume_mod += (volume == domains[domain_name]).astype(np.int) * self.sim.siteTypes[domain_name]

		# Set volume
		# Can be accelerated using serialization.
		for x in range(nx):
		    for y in range(ny):
		        for z in range(nz):
		            self.sim.lattice.setSiteType(x, y, z, int(volume_mod[x,y,z]))
		self.sim.hasBeenDiscretized = True

		# Register domain locations.
		## dict.fromkeys(domains, []) refers an identical empty list []. 
		self.locs   = dict.fromkeys(domains)
		for domain_name in domains:
		    self.locs[domain_name] = []

		for x in range(nx):
		    for y in range(ny):
		        for z in range(nz):
		            for domain_name in domains:
		                if (volume_mod[x,y,z] == self.sim.siteTypes[domain_name]):
		                    self.locs[domain_name].append((x,y,z))


		self.num_voxels = dict.fromkeys(domains)
		for domain_name in domains:
		    self.num_voxels[domain_name] = len(self.locs[domain_name])
                
		# unit converter
		NA           = 6.022e23
		experiment_volume_in_L = spacing_in_m * spacing_in_m * spacing_in_m * \
                                         1000 * nx * ny * nz
		self.per_uM = 1/(NA*1e-6*experiment_volume_in_L)


		# Extract the regions of surfaces and PSDs
		
		self.surf_voxel_locs = {}
		self.surf_voxel_prob = {}
		
		for k, v in surfaces.items():
			self.surf_voxel_locs[k] = np.nonzero(v > 0)
			self.surf_voxel_prob[k] = v[self.surf_voxel_locs[k][0], self.surf_voxel_locs[k][1], self.surf_voxel_locs[k][2] ]

		# print('self.memb_voxel_locs[0].shape[0]: ', self.memb_voxel_locs[0].shape[0])
		#print("len(locs['default'])   : ", len(self.locs['default']))
		#print("len(locs['cytoplasm']) : ", len(self.locs['cytoplasm']))
		#print("len(locs['psd'])       : ", len(self.locs['psd']))                           


                
	def number_per_1uM( self, domain_name ):
		num_voxels = self.num_voxels[domain_name]
		conc_in_uM = 1
		number_per_1uM = uM_to_num(conc_in_uM, num_voxels, self.spacing_in_m)
		return number_per_1uM

                
	def define_species(self, molecular_names):
		# Check arguments
		if isinstance(molecular_names, str):
		    molecular_names = [molecular_names]
		if not isinstance(molecular_names, (list, tuple)) :
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


	def reac_oneway_uM(self, reac, prod, rate, domain_name ):
		if isinstance(reac, tuple) and (len(reac) == 2):
		    kf  = self.per_uM * rate
		elif isinstance(reac, str):
		    kf  = rate
		else :
		    raise ValueError('Illigal values of reac in reac_oneway_uM.', reac)
		self.sim.modifyRegion(domain_name).addReaction(reactant=reac, product=prod, rate=kf)


	def reac_twoway_uM(self, reac, prod, rates, domain_name ):
		if isinstance(reac, tuple) and (len(reac) == 2):
		    kf  = self.per_uM * rates[0]
		elif isinstance(reac, str):
		    kf  = rates[0]
		else :
		    raise ValueError('Illigal values of reac in reac_twoway_uM.', reac)

		if isinstance(prod, tuple) and (len(prod) == 2):
		    kb  = self.per_uM * rates[1]
		elif isinstance(prod, str):
		    kb  = rates[1]
		else :
		    raise ValueError('Illigal values of prod in reac_twoway_uM.', prod)

		self.sim.modifyRegion(domain_name)\
			.addReaction(reactant=reac, product=prod, rate=kf)\
			.addReaction(reactant=prod, product=reac, rate=kb)


	def set_diffusion(self, molecular_name, rate, domain_name ):
		#
		self.sim.modifyRegion( domain_name ).setDiffusionRate(molecular_name, rate=rate)

	def add_molecule_uM(self, molecular_name, conc, domain_name):
		#
		number_per_1uM = self.number_per_1uM(domain_name)
		self.add_molecule( molecular_name, int(conc*number_per_1uM), domain_name )

	def add_molecule(self, molecular_name, number, domain_name):
		"""Add molecules to a domain.

		Args:
			molecular_name (str): Molecular name
			number (int): Number of molecules
			domain_name (str): Domain name

		Returns:
			(bool): True for success, False otherwise.
		"""
		particleNum=self.sim.particleMap[molecular_name]
		for i in random.sample(self.locs[domain_name], number):
		    self.sim.lattice.addParticle(i[0], i[1], i[2], int(particleNum))
		self.sim.customAddedParticleList.append((molecular_name, number))
		return True


	def add_surface_molecule(self, molecular_name, density, surface_name):
		"""Add surface molecules.

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


	def _check_arguments(self, volume, domains, spacing_in_m, surfaces):

		if not isinstance(volume, np.ndarray) or (volume.ndim != 3) or (volume.dtype not in NUMPY_INTEGERS):
			raise ValueError('volume must be a integer 3D np.ndarray.')
		elif not isinstance(domains, dict) :
			raise ValueError('domains must be dict.')
		elif not isinstance(surfaces, dict) :
			raise ValueError('surfaces must be dict.')
		elif not isinstance(spacing_in_m, float) :
			raise ValueError('spacing_in_m must be float.')
		else:
			return True

