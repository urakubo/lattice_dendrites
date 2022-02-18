.. include:: ../isonum.txt
.. include:: ../isogrk1.txt
.. |Ca2+| replace:: Ca\ :sup:`2+`

==============
Run simulation
==============

Models stored in LM-format files can be executed by Lattice Microbes (LM). LM is not a Python module, but is executed by the command prompt:

.. code-block:: bash

	$ mkdir results_photobleach
	$ cp models/photobleach.lm results_photobleach/0000.lm
	$ lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -f results_photobleach/0000.lm

The first command makes the directory 'results_photobleach', and the 
second command makes a copy of 'models/photobleach.lm' as 'results_photobleach/0000.lm', then it is executed by LM in the third command. This process can also be described using a Python script (41_single_run.py):


.. literalinclude:: ../../tutorial/1/41_single_run.py
   :language: python
   :linenos:
   :caption: 41_single_run.py

We can easily run simulation; however, in this simple run, we cannot handle any event functions such as photobleaching and synaptic input. The 'ConnectRun' class can handle those events by connecting multiple runs into a sequential run. Events are inserted just before the unit runs. In '42_connect_run.py', the class instance 'r' is created in Line 21, and 'r' accepts many variables such as a label volume file (label_volume_file; Line 22), template LM file (template_lm_file; Line 25), simulation result directory (output_dir; Line 26). The variable 'exec_periods' accepts the Python list that contains the simulation time of each unit run (Line 27), and 'exec_events' accepts the list of functions, each of which is called just before the unit run (Lines 28). 'event_params' sets parameters for the functions (Line 29). The class method 'exec' executes the defined connected run (Line 48).

.. literalinclude:: ../../tutorial/1/42_connect_run.py
   :language: python
   :linenos:
   :caption: 42_connect_run.py

The method 'exec' first calls 'null_event' (do nothing), and executes a 4-s pre-run. It then calls 'event_replace', which is described in Lines 6-18, and executes a 4-s post-run. The simulation results are stored in the directory 'output_dir' as '0000.lm' and '0001.lm'.

All event functions in the 'ConnectRun' class must accept three variables: lattice, sys_param, and event_param (Lines 6). 
The variable 'lattice' is a numpy 4D array that has a dimension of X |times| Y |times| Z |times| 16, the last of which denotes the slots of molecules (16) at each lattice site. The variable 'sys_param' conveys the system parameters such as the number of calls (sys_param['i']; Line 7), current simulation time (sys_param['time']; Line 8), a dict variable that specifies molecular species id (sys_param['species']; Line 9), and a label volume (sys_param['label volume']; Lines 14, 22). The variable 'event_param' contains the variable that is described in 'event_params' (Line 29). Here, we obtain the molecular ids of 'YFP' and 'bleached YFP' in the variables 'src' and 'dst', respectively. Then, the molecules 'YFP' are replaced with the molecules 'bleached YFP' if they are located in the labeled area (Lines 15-17).

All the event functions in the 'ConnectRun' must return the two variables: 'lattice' and 'event_param' (Lines 18). Both of them can be altered, and used for the subsequent simulation and event. In the case of |Ca2+| influx via NMDA receptors (Lines 32-37), the simulation time is decreased to 2 + 2 s, and 'event_replace' replaces the molecules 'inactive NMDAR' with the molecules 'active NMDAR'.


