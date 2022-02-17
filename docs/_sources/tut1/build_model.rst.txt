============
Build models
============
.. include:: ../isonum.txt
.. include:: ../isogrk1.txt
.. |Ca2+| replace:: Ca\ :sup:`2+`
.. |'''| replace:: \'\'\'
.. |um2| replace:: |mgr|\ m\ :sup:`2`

In this subsection, a LM model is built based on the voxelized shape of a spiny dendrite (31_build_model.py). 

1. Some HDF containers are loaded to define the cytosolic region (vol_cytosol; Line 21), PSD region (vol_PSD; Line 22), and the region of cell boundary (vol_bound; Line 23).

2. The dict variable 'domains' is set to designate the names of each integer in vol_cytosol (8-bit unsigned int; Line 25), and the dict variable 'surfaces' is set to designate the pairs of a name and voxelized surface areas (float; Line 27).

3. The BuildAnyShape class is called to register the domains and surfaces in the instance variable 'cell' (Line 32).

4. The instance variable 'cell' further incorporates the following properties regarding molecules: names, initial locations, diffusion, and interactions (set_molecules; Line 33). The contents of 'set_molecules' are described in the files 'set_molecule_FRAP.py' and 'set_molecule_Ca.py'. The Python language utilizes triple quotes (|'''|) to comment out a block of code. Thus users can comment out either of them (Lines 7-15).

5. Simulation timers are set under the instance 'sim' (Lines 33-37). Descriptions of setTimestep, setWriteInterval, setLatticeWriteInterval, and setSimulationTime are provided in `the instruction guide of LM <http://faculty.scs.illinois.edu/schulten/software_manuals/InstructionGuide.pdf>`_. Indeed, the BuildAnyShape class is a wrapper of a class of LM 'pyLM.RDME.RDMESimulation', all methods of which are preserved in the instance 'sim'. 

6. The overall setup is saved in the LM format file 'models/photobleach_yfp.lm' or 'models/Ca_influx.lm' (Lines 7-11 and 40-43).


.. literalinclude:: ../../tutorial/1/31_build_model.py
	:language: python
	:linenos:
	:caption: 31_build_model.py


|


.. image:: imgs/Scheme2.jpg
   :scale: 80%
   :align: center

|

We here recall two simulation targets: fluorescence recovery after photobleaching (FRAP; left in Figure above) and |Ca2+| influx via NMDA receptors (right in Figure above). To realize the FRAP simulation, the script below (set_molecule_FRAP.py) has a function 'set_molecules' to introduce two types of cytosolic molecules: 'YFP' and 'bleached YFP' (define_species; Lines 3-6), set their initial concentrations (add_molecule_uM; Lines 8-11), and diffuse them with a diffusion coefficient (set_diffusion; Lines 13-16). If we select to call the script 'set_molecule_FRAP.py' (Lines 7-15 in 31_build_model.py), these molecular settings are incorporated in the instance variable 'cell'.



.. literalinclude:: ../../tutorial/1/set_molecule_FRAP.py
	:language: python
	:linenos:
	:caption: set_molecule_FRAP.py


To simulate |Ca2+| influx via NMDA receptors, the script below (set_molecule_Ca.py) defines a variety of molecules and their properties. First, a cytosolic molecule '*Ca*' is created to represent |Ca2+| ions (Lines 4, 9). The initial number of '*Ca*' is set to be zero, which is declared using the number representation (add_molecule; Lines 13-15). Next, two PSD molecules '*active NMDAR*' and '*inactive NMDAR*' are created (Lines 5-6, 10), and 'inactive NMDAR' is set to have an initial number density of 300 /|um2| (add_surface_molecule; Lines 5-6, 16-17). Similary, the cell-boundary molecules '*Pump*' and '*Pump-Ca*' are created for |Ca2+| uptake (Lines 5-6, 10), and 'Pump' is distributed with a density of 100 /|um2| (add_surface_molecule; Lines 5-6, 16-17). Only Ca is set to be diffusible (Lines 19-21).


Then, molecular interactions are modeled (Figure below). The NMDA receptors are set to be rapidly activated at time 0 s, and this abrupt change will be introduced in the process of simulation. The '*active NMDAR*' mediates |Ca2+| influx. This is formalized by a first-order reaction, with a rate constant of '*k_channel_Ca*' (reac_oneway_uM; Lines 30, 33-34). The '*active NMDAR*' is gradually converted to its inactivated form '*inactive NMDAR*' (first-order rate constant: '*k_nmdar_deact*'; Lines 31, 33-34). Cytosolic '*Ca*' is uptaken by a pump to the extracellular space (Figure below). In this process, '*Ca*' binds to '*Pump*' in a reversible manner (reac_twoway_uM; Lines 36-39), and the |Ca2+|-bound pump '*Pump-Ca*' releases the |Ca2+| to the extracellular space (first-order rate constant: '*k_pump_Ca*'; Lines 32-34). These interactions summarize a life cycle of spine |Ca2+| that is a cause of synaptic plasticity.

|

.. image:: imgs/Scheme3.jpg
   :scale: 80%
   :align: center

.. literalinclude:: ../../tutorial/1/set_molecule_Ca.py
	:language: python
	:linenos:
	:caption: set_molecule_Ca.py

|

The successful distribution of molecules can be confirmed by executing '32_show_molecule.py' (Figure below). You can find this script in $LD_DIRECTORY/tutorial/1 .

.. image:: imgs/Initials.jpg
   :scale: 100%
   :align: center

