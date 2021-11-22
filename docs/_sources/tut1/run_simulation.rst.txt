==============
Run simulation
==============

We would like to keep "models/ball_and_stick.lm" as a template, and store the simulation reults in "results/run_ball_and_stick.lm". To do this, we first copy the file then execute lm simulation through the command prompt:

.. code-block:: bash

	$ mkdir results
	$ cp models/ball_and_stick.lm results/run_ball_and_stick.lm
	$ lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -f results/run_ball_and_stick.lm


We can also write a Python script to realize the above:

.. literalinclude:: ../../tutorial/1/41_single_run.py
   :language: python
   :linenos:
   :caption: 41_single_run.py


We can easily run such a single simulation; however, it would not satisfy an aim of simulation.
For example, one may want to see a fluorescence recovery after photobleaching (FRAP).
In this case, a preparatory run should be followed by an instantaneous action, i.e., photobleach.
Then, the recovery process would be observed.
This is raelized by the following script.

.. literalinclude:: ../../tutorial/1/42_connect_run.py
   :language: python
   :linenos:
   :caption: 42_connect_run.py


That is all for simulation run.
