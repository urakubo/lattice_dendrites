==============
Run simulation
==============

We would like to keep "models/ball_and_stick.lm" as a template, and store the simulation reults in "results/run_ball_and_stick.lm". To do this, we first copy the lm file then execute the simulation through the command prompt:

.. code-block:: bash

	$ mkdir results
	$ cp models/ball_and_stick.lm results/run_ball_and_stick.lm
	$ lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -f results/run_ball_and_stick.lm


We can also write a Python script like:

.. literalinclude:: ../../tutorial/1/41_single_run.py
   :language: python
   :linenos:
   :caption: 41_single_run.py


We can easily run such a simulation; however, it would not satisfy aims of simulation.
For example, one may want to see a fluorescence recovery after photobleaching (FRAP).
In this case, a preparatory run should be followed by an instantaneous action, i.e., photobleach.
Then, the recovery process should be simulated.
To handle such events, we can use the "ConnectRun" class as follows:

.. literalinclude:: ../../tutorial/1/42_connect_run.py
   :language: python
   :linenos:
   :caption: 42_connect_run.py


That is all for simulation run.
