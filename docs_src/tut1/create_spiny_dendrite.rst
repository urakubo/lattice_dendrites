.. include:: ../isonum.txt
.. include:: ../isogrk1.txt

==============
Create a shape
==============

1. We would first draw a blueprint of a spiny dendrite as a combination of geometric shapes (Figure below). A sphere with a radius of 0.25 |mgr| m represents a spine head, part of which is labeled as a postsynaptic density (PSD; red colored area). The spine head has a cylindrical spike neck (radius: 0.1 |mgr| m, length 1.0 |mgr| m), which is connected to a parent dendrite that also has a cylindrical shape (radius: 0.5 |mgr| m, length 1.8 |mgr| m). Not only the contour of spiny dendrite, we also introduce two types of intracellular organelles. Endoplasmic reticulum (ER) is set as thin cylinders (radius: 0.08 |mgr| m, length 1.8 |mgr| m; blue colored objects), and a mitochondrion is set as a thich cylinder (radius: 0.2 |mgr| m, length 1.8 |mgr| m; yellow colored objects). 

|

.. image:: imgs/Scheme.jpg
   :scale: 60%
   :align: center


|

2. The designed shape is embedded in a voxel space for simulation (Lines 1-35, 11_create_dend.py). First, a voxel space of 96 |times| 60 |times| 96 voxels (20 nm/voxel) is set as a numpy 3D array 'vol_dend' (Line 14). Then, the spine head is created using a Python module function morphology.ball (skimage; Line 10; 0.25 |mgr| m, 12 voxels), which is added to 'vol_dend' (add_shape; Line 15). Similarly, we make the spine neck and dendrite as cylinders (Lines 11 and 12, respectively), and add them to 'vol_dend' (Lines 16 and 17, respectively). The in-house function 'add_shape,' simply overlays overwrapped regions, and represents the filled areas as 1 (cytosolic region) and the void areas as 0 (extracellular region). Thus, we have already built the contour of the spiny dendrite in 'vol_dend.' Similarly, the PSD, Mitochondrion, and ER are embedded in vol_psd, vol_mito, and vol_er, respectively.

3. Here, the volume sizes need to be corrected. LM can only simulate a volume with a size of a multiple of 32 |times| 32 |times| 32 voxels, whereas the current voxel size is 96 |times| 60 |times| 96 voxels. We thus execute a utility function 'lmpad' (Lines 39-42). The voxel sizes are automatically expanded to the multiple of 32 (96 |times| 64 |times| 96 voxels in this case). Users can of course set a multiple of 32 voxels from the beginning.



.. literalinclude:: ../../tutorial/1/11_create_dend.py
   :language: python
   :linenos:
   :caption: 11_create_dend.py


4. The surface of each object is calculated.  

.. literalinclude:: ../../tutorial/1/12_show_dend.py
   :language: python
   :linenos:
   :caption: 12_show_dend.py

.. image:: imgs/ball_and_stick.png
   :scale: 50%
   :align: center

That is all.

