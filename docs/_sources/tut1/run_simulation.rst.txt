==============
Run simulation
==============

Models stored in LM-format files can be executed by Lattice Microbes (LM). LM is not a Python module, but is executed by the command prompt:

.. code-block:: bash

	$ mkdir results_photobleach
	$ cp models/photobleach.lm results_photobleach/0000.lm
	$ lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -f results_photobleach/0000.lm

The above commands make a copy of 'models/photobleach.lm' as 'results_photobleach/0000.lm', and it is executed by LM. This process can also be described using a Python script, as follow (41_single_run.py):

.. literalinclude:: ../../tutorial/1/41_single_run.py
   :language: python
   :linenos:
   :caption: 41_single_run.py

We can easily run simulation; however, in this simple run, we cannot handle any event functions, such as photobleaching and synaptic input. The "ConnectRun" class enables them by connecting multiple runs into a sequential run. Events can be inserted in between the unit runs. The ‘42_connect_run.py’ is an example:

.. literalinclude:: ../../tutorial/1/42_connect_run.py
   :language: python
   :linenos:
   :caption: 42_connect_run.py


That is all for simulation run.
