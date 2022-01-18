# 
# University of Illinois Open Source License
# Copyright 2008-2016 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Michael J. Hallock and Joseph R. Peterson
#
# 

import lm, math, os
from .LMLogger import *
from .ipyInterface import *
try:
	from tqdm import tqdm
except:
	def tqdm(x,ascii=False):
		return x


## @class CMESimulation
# A CME simulation that contains reactions and species
class CMESimulation:
	## A constructor fo the CMESimulation
	# @param volume The reaction vessel volume (in Litres).  Specifying 'None' signifies that the user has already accounted for volume in rate constants.
	# @param name The name of the CME simulation; Default: "unnamed"
	# @return self
	def __init__(self,volume = None, name="unnamed"):

		self.particleMap={}
		self.species_id=[]
		self.initial_counts={}
		self.reactions=[]
		self.parameters={}
		self.volume = volume
		self.name = name

		self.replicates = []   # Replicates that have been run
		self.filename   = None # The filename used when running

	## Define a specie/s of particles that exist in the simulation
	# @param self
	# @param species A list of species to add to the simulation
	def defineSpecies(self, species):
		for s in species:
			self.species_id.append(s)
			self.particleMap[s]=len(self.species_id)
			self.initial_counts[s]=0

	## Add a specified number of particles of the specified type to the  specified region
	# @param self
	# @param species The name of the specie to add particles to
	# @param count The number of that particle to start with (default 1)
	def addParticles(self, species='unknown', count=1):
		try:
			particleNum=self.particleMap[species]
		except KeyError:
			particleNum=1
			LMLogger.warn('In CME.addParticles, couldn\'t find particle of type "%s" in map (is it previously defined with \'CME.defineSpecies(...)\'?).  Shouldn\'t happen.',species)
		self.initial_counts[species]+=count

	## Add a concentration of particles of the specified type to the simulation
	# @param self
	# @param species The name of the specie to add
	# @param concentration The concentration of the species (Molar).  Particle count is rounded to nearest integer.
	# @throws An error if the volume is not specified
	def addConcentration(self, species='unknown', conc=0.0):
		if self.volume == None:
			LMLogger.error('In CME.addConcentration, a volume must be specified in CME constructor.')
			raise Exception("No volume specified.")
		try:
			particleNum=self.particleMap[species]
		except KeyError:
			particleNum=1
			LMLogger.warn('In CME.addConcentration, couldn\'t find particle of type "%s" in map (is it previously defined with \'CME.defineSpecies(...)\'?).  Shouldn\'t happen.',species)
		self.initial_counts[species]= int(round(conc*self.volume*6.022e23))


	## Adds a 0th, 1st or 2nd order reaction
	# @param self 
	# @param reactant The list or reactants
	# @param product The list of products
	# @param rate The rate of reaction
	def addReaction(self, reactant, product, rate):
		if rate <= 0.0:
			LMLogger.error("In CME.addReaction, the rate must be a positive number")
		if (reactant,product,rate) in self.reactions:
			LMLogger.warning("Reaction already in model: %s -> %s : %f"%(reactant,product,rate))
			return
		self.reactions.append((reactant,product,rate))

	## Return the Lattice Microbes ReactionModel object for fine-tuning
	# @param self
	# @return The reaction model for the simulation
	def buildReactionModel(self):
		rm=lm.ReactionModel()

		numReactions=len(self.reactions)
		LMLogger.info("number of reactions = %d", numReactions)
		rm.set_number_reactions(numReactions)

		numSpecies=len(self.species_id)
		LMLogger.info("number of species = %d", numSpecies)
		rm.set_number_species(numSpecies)

		for s in self.species_id:
			LMLogger.debug("\t set initial counts of %s (id %d) to %d", s, self.particleMap[s], self.initial_counts[s])
			rm.add_initial_species_count(self.initial_counts[s])

		# Build reaction matricies
		rxtypes = [0] * numReactions;
		rxconst = [0] * numReactions;
		stoich  = [0] * numSpecies * numReactions;
		depmat  = [0] * numSpecies * numReactions;

		rnum=0
		for rx in self.reactions:
			reactant=rx[0]
			product=rx[1]
			rate=rx[2]
			# Scale rates by volume, if specified
			if self.volume != None:
				if isinstance(reactant, tuple):
					# Second order	
					if reactant[0] != reactant[1]:
						rate /= (6.022e23*float(self.volume))
					elif reactant[0] == reactant[1]:
						# Second order self-reaction
						rate /= (6.022e23*float(self.volume))
				else:
					# Zeroth order
					if reactant == '':
						rate *= 6.022e23*float(self.volume)
					# First order
					# do nothing...


			LMLogger.debug("\t rx %d is %s -> %s at rate %g", rnum, reactant, product, rate)
			rxconst[rnum]=rate
			if(isinstance(reactant, tuple)):
				# Second Order
				if reactant[0] != reactant[1]:
					rxtypes[rnum]=2
					for sr in reactant:
						depmat[rnum + numReactions*(self.particleMap[sr]-1)]=1
						stoich[rnum + numReactions*(self.particleMap[sr]-1)] -= 1
				elif reactant[0] == reactant[1]:
				# Second Order Self Reaction
					rxtypes[rnum]=3
					depmat[rnum + numReactions*(self.particleMap[reactant[0]]-1)]=1
					stoich[rnum + numReactions*(self.particleMap[reactant[0]]-1)] -= 2
			else:
				# Zeroth order
				if reactant == '':
					rxtypes[rnum] = 0
					#stoich[rnum + numReactions*(self.particleMap[reactant]-1)] += 0
				else:
				# First order
					rxtypes[rnum]=1
					depmat[rnum + numReactions*(self.particleMap[reactant]-1)]=1
					stoich[rnum + numReactions*(self.particleMap[reactant]-1)] -= 1

			if(isinstance(product, tuple)):
				for sp in product:
					stoich[rnum + numReactions*(self.particleMap[sp]-1)] += 1
			else:
				if product != '':
					stoich[rnum + numReactions*(self.particleMap[product]-1)] += 1

			rnum += 1


		LMLogger.debug("rxtypes: %s", rxtypes)
		LMLogger.debug("rxtypes: %s", rxconst)

		# Populate reaction object data
		for r in range(numReactions):
			rm.add_reaction()
			rm.mutable_reaction(r).set_type(rxtypes[r])
			rm.mutable_reaction(r).add_rate_constant(rxconst[r])

		LMLogger.debug("stoich:  %s", stoich)
		for x in range(len(stoich)):
			rm.add_stoichiometric_matrix(stoich[x])

		LMLogger.debug("depmat:  %s", depmat)
		for x in range(len(depmat)):
			rm.add_dependency_matrix(depmat[x])

		return rm

	# Set a known random seed
	# @param self
	# @param seed The seed value
	def setRandomSeed(self, seed):
		self.parameters['seed']=str(seed)


	## Set the simulation state write-to-disk interval	
	# @param self
	# @param time Time length between writes
	def setWriteInterval(self, time):
		self.parameters['writeInterval']=str(time)

	## Set the total simulation time
	# @param self
	# @param time Time length of the simulation
	def setSimulationTime(self, time):
		self.parameters['maxTime']=str(time)

	# set the simulation time step
	# @param self
	# @param time The length of simulation timestep for CME
	def setTimestep(self, time):
		self.parameters['timestep']=str(time)

	# ## Check that each simulation parameter is reasonable
	# # @param self
	# def check(self):
	# 	pass

	## Create an HDF5 version of the simulation amenable for later running or stand-alone running
	# @param self
	# @param filename The filename to save the simulation setup in
	def save(self, filename):
		lm.SimulationFile.create(filename)
		f=lm.SimulationFile(filename)

		rm=self.buildReactionModel()
		f.setReactionModel(rm)

		for key in self.parameters:
			LMLogger.debug("Setting parameter %s = %s", key, self.parameters[key])
			f.setParameter(key, self.parameters[key])
		LMLogger.debug("SpeciesNames: %s", ",".join(self.species_id))
		f.setParameter("speciesNames", ",".join(self.species_id))

		f.close()

	## Run the simulation using the specified solver the specified amount of time
	# @param self
	# @param filename The HDF file to write to
	# @param method The class name for the solver to use (e.g., lm::cme::GillespieDSolver")
	# @param replicates The number of replicates to serially run
	# @param seed A seed for the random number generator to use when running the simulation; None indicates default
	def run(self, filename, method, replicates=1, seed=None, cudaDevices=None, checkpointInterval=0):

		if cudaDevices is None:
			cudaDevices = [0]

		if seed is not None:
			f = lm.SimulationFile(filename)
			f.setParameter("seed",str(seed))
			f.close()
		for r in tqdm(range(1, replicates+1),ascii=True):
			lm.runSimulation(filename, r, method, cudaDevices, checkpointInterval)
			# update internal state with replicates that have been run
			self.replicates.append(r)
		# Update the filename
		self.filename = filename
	
	## Run the simulation using a call to mpirun with the given options
	# @param self
	# @param filename The HDF file to write to
	# @param method The class name for the solver to use (e.g., lm::cme::GillespieDSolver")
	# @param replicates The number of replicates to serially run
	# @param driver The program to execute the parallel run, e.g. "mpirun", "aprun", "ibrun", etc.
	# @param ppe The number of processing elements to use 
	# @param seed A seed for the random number generator to use when running the simulation; None indicates default
	def runMPI(self, filename, method, replicates=1, driver="mpirun", ppe=1, seed=None):
		if seed is not None:
			f = lm.SimulationFile(filename)
			f.setParameter("seed",str(seed))
			f.close()
		# Write out system command
		repStr = "1"
		if replicates != 1:
			repStr = "1-%d"%replicates
		cmdStr = '%s -np %d  lm -ws -r %s -sl "%s"  -f %s'%(driver, ppe, repStr, method, filename)
		LMLogger.debug("Running MPI LM job with command:\n\t%s"%(cmdStr))

		# Actually run command and check return code
		status = os.system(cmdStr)
		
		# Check that LM executed correctly
		if status != 0:
			LMLogger.error("Failed running MPI (status: %d) job running command:\n\t%s\n"%(status, cmdStr))

		# Update the internal state for replicates that have been run
		for r in range(1,replicates):
			self.replicates.append(r)
		# Update the filename
		self.filename = filename

	# Run a simulation with a given solver
	# @param self
	# @param filename The HDF file to write to
	# @param solver An MESolver object
	# @param replicates The number of replicates to serially run
	def runSolver(self, filename, solver, replicates=1, cudaDevices=None, checkpointInterval=0):
		if cudaDevices is None:
			cudaDevices = [0]
		LMLogger.debug("Running Custom CMESolver: %s"%(solver))
		for r in tqdm(range(1, replicates+1),ascii=True):
			LMLogger.debug("  Running replicate: %d"%(r))
			lm.runSolver(filename, r, solver, cudaDevices, checkpointInterval)
			# update internal state with replicates that have been run
			self.replicates.append(r)
		# Update the filename
		self.filename = filename


	# Create a representation that is loadable in iPython
	# @return A string containing an HTML object that can be displayed in the Jupyter notebook
	def _repr_html_(self):
		# Write out major parameters
		s  = ""
		s += "<h1> CME Simulation: %s</h1>"%(self.name)
		if self.volume or len(self.parameters) > 0:
			s += "<h2> Simulation Parameters </h2>"
			if self.volume:
				s += "Volume: %f<br/>"%float(self.volume)
			for k,v in self.parameters.items():
				s += "%s: %s</br>"%(str(k), str(v))
		# Write completed simulations
		if len(self.replicates) > 0:
			s += "<h2> Simulations </h2>"
			s += "Filename: " + str(self.filename) + "<br/>"
			s += "Completed replicates: " + str(self.replicates) + "<br/>"
		# Write out species
		if len(self.species_id) > 0:
			s += "<h2> Species </h2>"
			s += writeTable(["Specie","Particle ID","Initial Count"], [(specie, self.particleMap[specie], self.initial_counts[specie]) for specie in self.species_id])		
		# Write out reactions
		if len(self.reactions) > 0:
			s += "<h2> Reaction Model </h2>"
			s += writeTable(["Reaction","Rate","Units"], [getReactionString(rxn[0],rxn[1],rxn[2]) for rxn in self.reactions])
		# HTML string for representing the model
		return s

