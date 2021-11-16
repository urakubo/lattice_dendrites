=======================
Visualization 2: mayavi
=======================

Some users need to further analyze simulation results in a spatio-temporal manner, e.g., for 3D visualization. The class ConnectAnalysis allows to use user-defined function for a time series of 4D images (3D space + 16 slots for molecules).

.. literalinclude:: ../../tutorial/1/61_video_photobleach.py
   :language: python
   :linenos:
   :caption: 61_video_photobleach.py

ConnectAnalysis calls the function "event" at each time of recording, and it produces a mayavi-based video as follows:

.. video:: _static/photobleach.mp4

That is all for visualzation 2.
