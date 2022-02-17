==============
Run simulation
==============

Models stored in LM-format files can be executed by Lattice Microbes (LM). LM is not a Python module, but a standalone software that can be executed through the command prompt:

.. code-block:: bash

	$ mkdir results_photobleach
	$ cp models/photobleach.lm results_photobleach/0000.lm
	$ lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -f results_photobleach/0000.lm

These commands make a copy of 'models/photobleach.lm' as 'results_photobleach/0000.lm', which is executed by lm. This process can also be described using a Python script.

.. literalinclude:: ../../tutorial/1/41_single_run.py
   :language: python
   :linenos:
   :caption: 41_single_run.py


We can easily run simulation; however, using this simple run, we cannot introduce any events function, such as photobleaching and synaptic input. To handle such events, the "ConnectRun" class enables the is provided as used in '42_connect_run.py'.

.. literalinclude:: ../../tutorial/1/42_connect_run.py
   :language: python
   :linenos:
   :caption: 42_connect_run.py


That is all for simulation run.
