.. include:: ../isonum.txt
.. include:: ../isogrk1.txt

==============
Create a shape
==============

|

.. image:: imgs/Scheme.jpg
   :scale: 60%
   :align: center


|


We would first draw a blueprint of a schematics shape of spiny dendrite (Figure above), and the designed shape is embedded in a voxel space.

1. a voxel space of 96 |times| 60 |times| 96 voxels (20 nm/voxel) is set as a numpy 3D array 'vol_dend' (Line 14 in 11_create_dend.py). Then, we use a Python module function morphology.ball (skimage) to create a spine head (Line 10) as a sphere with a radius of 0.25 |mgr| m (12 voxels). This sphere is added to the voxel space 'vol_dend' (add_shape; Line 15). Similarly, we make spine neck and dendrite as cylinders (Lines 11 and 12, respectively), and add them to 'vol_dend' (add_shape; Lines 16 and 17, respectively). The in-house function 'add_shape,' simply overlays overwrapped regions, and represents the filled areas as 1 (cytosolic region) and the void areas as 0 (extracellular region). Thus, we have already built the contour of the spiny dendrite in 'vol_dend.'



.. literalinclude:: ../../tutorial/1/11_create_dend.py
   :language: python
   :linenos:
   :caption: 11_create_dend.py

.. literalinclude:: ../../tutorial/1/12_show_dend.py
   :language: python
   :linenos:
   :caption: 12_show_dend.py

.. image:: imgs/ball_and_stick.png
   :scale: 50%
   :align: center

That is all.

