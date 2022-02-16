=============
Set reactions
=============


In this subsection, a LM model is built based on the voxelized shape of a spiny dendrite (31_build_model.py). 

1. We first load the HDF containers that has the cytosolic region (vol_cytosol; Line 18), PSD region (vol_PSD; Line 19), and the region of cell boundary (vol_bound; Line 20).

2. We set the dict variable 'domains' to designate the names of each integer in vol_cytosol (8-bit unsigned int; Line 23), and the dict variable 'surfaces' to designate the pairs of a name and voxelized surface areas (64-bit float; Line 24).

3. We then call the BuildAnyShape class (Line 29) to register the domains and surfaces.

4. We further set molecular properties, such as names, initial locations, diffusion, and interactions (set_molecules; Line 30). They are described in the files '32_FRAP.py' and '33_CaSignal.py', either of which can be incorporated (Lines 7-11).

5. We set simulation timers (Lines 33-37), and save the LM model (Lines 40-43).



.. literalinclude:: ../../tutorial/1/31_build_model.py
	:language: python
	:linenos:
	:caption: 31_build_model.py


|

We here recall the two simulation targets: fluorescence recovery after photobleaching (FRAP; left in Figure below) and Ca2+ influx through NMDA receptors (right in Figure below). 

|



.. image:: imgs/Scheme2.jpg
   :scale: 80%
   :align: center

|

Two parts can be replace as follows to exmine the reactions ('32_FRAP.py'):

.. literalinclude:: ../../tutorial/1/32_FRAP.py
	:language: python
	:linenos:
	:caption: 32_FRAP.py

|

33_Ca.py

|

.. literalinclude:: ../../tutorial/1/33_CaSignal.py
	:language: python
	:linenos:
	:caption: 33_CaSignal.py


That is all for simulation setup.
