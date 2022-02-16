============
Build models
============

.. |Ca2+| replace:: Ca\ :sup:`2+`

In this subsection, a LM model is built based on the voxelized shape of a spiny dendrite (31_build_model.py). 

1. Some HDF containers are loaded to define the cytosolic region (vol_cytosol; Line 18), PSD region (vol_PSD; Line 19), and the region of cell boundary (vol_bound; Line 20).

2. The dict variable 'domains' is set to designate the names of each integer in vol_cytosol (8-bit unsigned int; Line 23), and the dict variable 'surfaces' is set to designate the pairs of a name and voxelized surface areas (64-bit float; Line 24).

3. The BuildAnyShape class is called to register the domains and surfaces in the instance variable 'cell' (Line 29).

4. The instance variable 'cell' further incorporates the following properties regarding molecules: names, initial locations, diffusion, and interactions (set_molecules; Line 30). The entities are described in the files '32_FRAP.py' and '33_CaSignal.py', either of which can be incorporated (Lines 7-11).

5. Simulation timers are set under the instance 'sim' (Lines 33-37). Descriptions of setTimestep, setWriteInterval, setLatticeWriteInterval, and setSimulationTime are provided in `the instruction guide of LM <http://faculty.scs.illinois.edu/schulten/software_manuals/InstructionGuide.pdf>`_. Indeed, the BuildAnyShape class is a wrapper of a class of LM 'pyLM.RDME.RDMESimulation', all methods of which are preserved in the instance 'sim'. 

6. The overall setup is saved in the file 'models/photobleach_yfp.lm' or 'models/Ca_influx.lm' (Lines 7-11 and 40-43).


.. literalinclude:: ../../tutorial/1/31_build_model.py
	:language: python
	:linenos:
	:caption: 31_build_model.py


|


.. image:: imgs/Scheme2.jpg
   :scale: 80%
   :align: center

|

We here recall two simulation targets: fluorescence recovery after photobleaching (FRAP; left in Figure above) and |Ca2+| influx through NMDA receptors (right in Figure above). To realize the FRAP simulation, the script below (32_FRAP.py) has a function 'set_molecules' to introduce two types of cytosolic molecules: 'YFP' and 'bleached YFP' (define_species; Lines 3-6), set their initial concentrations (add_molecule_uM; Lines 8-11), and diffuse them with a diffusion coefficient (set_diffusion; Lines 13-16). If we select to call the script '32_FRAP.py' (Lines 7-11 in 31_build_model.py), the definitions of the molecules are incorporated in the instance variable 'cell'.



.. literalinclude:: ../../tutorial/1/set_molecule_FRAP.py
	:language: python
	:linenos:
	:caption: set_molecule_FRAP.py

|

In the case of |Ca2+| influx via NMDA receptors, the script '33_CaSignal.py' defines a variety of molecules and their properties. A cytosolic molecule 'Ca' is created to represent |Ca2+| ions. ‘Pump’ and ‘Pump-Ca’ are created for the uptake of |Ca2+| ions at the cell boundary.

|

.. literalinclude:: ../../tutorial/1/set_molecule_Ca.py
	:language: python
	:linenos:
	:caption: set_molecule_Ca.py


That is all for simulation setup.
