.. include:: ../isonum.txt
.. include:: ../isogrk1.txt
.. |Ca2+| replace:: Ca\ :sup:`2+`

======================
Visualization 2: video
======================
Some users may want to further analyze simulation results in a spatio-temporal manner, e.g., 3D visualization. The script '61_make_video.py' introduces a method to make videos of molecular dynamics using the visualization software Mayavi and the video generation software ffmpeg. In this script, the ConnectAnalysis class offers the use of user-defined function to access all dimension of the simulation data, i.e., a time series of 4D images (3D space + 16 molecular slots; Lines 55-60). The time-developing 4D images are repeatedly passed to the function 'event' (Lines 34-52).

The function 'event' has three variables: the lattice, sys_param, and usr_param. In 'event', the function 'get_loc_molecules' provides the locations of a specified molecule in the specified lattice (Lines 10, 38, 40), and the molecular locations are plotted in the pre-defined 3D space (Lines 31, 41-43). The visualized molecular distribution is saved as XXXX.png in 'output_dir' (Lines 24, 47-48). This process is repeated in a time-development manner.    

.. literalinclude:: ../../tutorial/1/61_make_video.py
   :language: python
   :linenos:
   :caption: 61_make_video.py

|

Finally, the video generation software 'ffmpeg' assembles the sequential images into a single video (Lines 63-69). The generated videos are shown below. In the case of FRAP (Lines 9-13), we can observe the re-filling of YFP (yellow points) to the target spine after the photobleaching.

.. video:: _static/YFP.mp4

|

In the case of |Ca2+| influx via NMDA receptors (Lines 16-20), we can observe the spread of |Ca2+| ions (blue points) from the stimulated spine.

.. video:: _static/Ca.mp4

