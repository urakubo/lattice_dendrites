from __future__ import print_function
from __future__ import division
from pyLM import *
import numpy as np
import os
import random


class buildAnyShape:
    def __init__(self, sim, volume, domains, voxel_membrane_area, voxel_PSD, voxel_ER_area):

		"""Exceptions are documented in the same way as classes.

		The __init__ method may be documented in either the class level
		docstring, or as a docstring on the __init__ method itself.
		Either form is acceptable, but the two should not be mixed.
		Note:
		    Do not include the `self` parameter in the ``Args`` section.

		Args:
		    msg (str): Human readable string describing the exception.
		    code (:obj:`int`, optional): Error code.

		Attributes:
		    msg (str): Human readable string describing the exception.
		    code (int): Exception error code.

		"""

        """
        
		sim.siteTypes['domain name'] may be different from volume_id.
		
		param sim           RDMESimulation object
		param volume        Three-dimentional array that indicates volume_id
		param domains_ dict {'domain name', volume_id}
		param voxel_membrane_area  Three-dimentional array that specifies membrane_voxels > 0
		param voxel_PSD   Three-dimentional array that specifies PSD_voxels > 0

		return self
		"""

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
                    sim.lattice.setSiteType(x, y, z, volume_mod[x,y,z])
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
        self.memb_voxel_locs = np.nonzero(voxel_membrane_area > 0)
        self.memb_voxel_prob = voxel_membrane_area[self.memb_voxel_locs[0],\
                                             self.memb_voxel_locs[1],\
                                             self.memb_voxel_locs[2] ]
        self.PSD_voxel_locs = np.nonzero(voxel_membrane_area*voxel_PSD > 0)
        self.PSD_voxel_prob = voxel_membrane_area[self.PSD_voxel_locs[0],\
                                             self.PSD_voxel_locs[1],\
                                             self.PSD_voxel_locs[2] ]
        self.ER_voxel_locs = np.nonzero((volume == domains['cytoplasm'])*(voxel_ER_area > 0))
        self.ER_voxel_prob = voxel_ER_area[self.ER_voxel_locs[0],\
                                             self.ER_voxel_locs[1],\
                                             self.ER_voxel_locs[2] ]
        # print('self.memb_voxel_locs[0].shape[0]: ', self.memb_voxel_locs[0].shape[0])

        
        #print("len(locs['default'])   : ", len(self.locs['default']))
        #print("len(locs['cytoplasm']) : ", len(self.locs['cytoplasm']))
        #print("len(locs['psd'])       : ", len(self.locs['psd']))                           

    def addCytosolicMolecules(self, molecular_name, molecular_number, domain_name):
        ## 
        ## Distribute specified molecules randomly in specified domains
        ##
        #
        # @param molecular_name
        # @param molecular_number
        # @param domain_name
        # @return self

        particleNum=self.sim.particleMap[molecular_name]
        for i in random.sample(self.locs[domain_name], molecular_number):
            self.sim.lattice.addParticle(i[0], i[1], i[2], particleNum)
        self.sim.customAddedParticleList.append((molecular_name, molecular_number))

    def addMembraneMolecules(self, molecular_name, density):
        ##
        particleNum=self.sim.particleMap[molecular_name]
        molecular_numbers = np.random.binomial(density, self.memb_voxel_prob)
        # print('np.sum(molecular_numbers): ', np.sum(molecular_numbers))
        for x, y, z, num in zip(self.memb_voxel_locs[0], \
                                self.memb_voxel_locs[1], \
                                self.memb_voxel_locs[2], \
                                molecular_numbers):
            # print('x,y,z, num: ', x,y,z, num)
            for i in range(num):
                self.sim.lattice.addParticle(x, y, z, particleNum)

    def addPSDMolecules(self, molecular_name, density):
        ##
        particleNum=self.sim.particleMap[molecular_name]
        molecular_numbers = np.random.binomial(density, self.PSD_voxel_prob)
        # print('np.sum(molecular_numbers): ', np.sum(molecular_numbers))
        for x, y, z, num in zip(self.PSD_voxel_locs[0], \
                                self.PSD_voxel_locs[1], \
                                self.PSD_voxel_locs[2], \
                                molecular_numbers):
            # print('x,y,z, num: ', x,y,z, num)
            for i in range(num):
                self.sim.lattice.addParticle(x, y, z, particleNum)


    def addERMolecules(self, molecular_name, density):
        ##
        particleNum=self.sim.particleMap[molecular_name]
        molecular_numbers = np.random.binomial(density, self.ER_voxel_prob)
        # print('np.sum(molecular_numbers): ', np.sum(molecular_numbers))
        for x, y, z, num in zip(self.ER_voxel_locs[0], \
                                self.ER_voxel_locs[1], \
                                self.ER_voxel_locs[2], \
                                molecular_numbers):
            # print('x,y,z, num: ', x,y,z, num)
            for i in range(num):
                self.sim.lattice.addParticle(x, y, z, particleNum)

