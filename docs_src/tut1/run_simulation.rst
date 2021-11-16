==============
Run simulation
==============

We would lienk to keep "models/ball_and_stick.lm" as a template, and store the simulation reults in "results/run_ball_and_stick.lm". We can realize it from the command prompt as follows.

.. code-block:: bash

	$ mkdir results
	$ cp models/ball_and_stick.lm results/run_ball_and_stick.lm
	$ lm -r 1 -sp -sl lm::rdme::MpdRdmeSolver -f results/run_ball_and_stick.lm


We can also arrage the above as a Python script as follows.

.. literalinclude:: ../../tutorial/1/41_single_run.py
   :language: python
   :linenos:
   :caption: 41_single_run.py


In many cases, such a single run would not satisfy the aim of simulation. The repeat of runs is wrapped in the following code.


.. literalinclude:: ../../tutorial/1/42_repeat_run.py
   :language: python
   :linenos:
   :caption: 42_repeat_run.py


That is all for simulation run.
