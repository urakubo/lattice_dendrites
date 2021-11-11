==============
Run simulation
==============

The following code 

.. code-block:: python
	:linenos:

	import h5py
	import numpy as np
	import os, sys, shutil
	import subprocess as s

	filename_lm     = 'morph/ball_and_stick.lm'
	filename_prerun = 'results/prerun_ball_and_stick.lm'

	if os.path.isfile(filename_prerun):
	    print('Prerun file was removed.')
	    os.remove(filename_prerun)

	print('Create prerun file: ', filename_prerun)
	shutil.copy(filename_lm, filename_prerun)

	com = ['lm','-r', '1', '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_prerun]
	print(' '.join(com))
	s.call(com)


.. code-block:: python
	:linenos:

	gpu_id = '0,1'
	num_GPUs = '2'
	num_CPUs = '2'
	com = ['lm','-g', gpu_id, '-gr', num_GPUs,'-cr', num_CPUs, '-sl','lm::rdme::MGPUMpdRdmeSolver','-f', filename_prerun]

今回は、index.rstの中身は必要ないので、一旦全て削除します。
