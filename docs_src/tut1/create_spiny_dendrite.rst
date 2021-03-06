.. include:: ../isonum.txt
.. include:: ../isogrk1.txt

==============
Create a shape
==============

1. We first draw a blueprint to make a spiny dendrite as a combination of geometric shapes (Figure below). A sphere (radius: 0.25 |mgr|\ m) represents a spine head, which is partially labeled as a postsynaptic density (PSD; red colored area). The spine head has a cylindrical spike neck (radius: 0.1 |mgr|\ m, length 1.0 |mgr|\ m), and it is further connected to a parent dendrite that also has a cylindrical shape (radius: 0.5 |mgr|\ m, length 1.8 |mgr|\ m). Not only the contour of spiny dendrite, we also introduce two types of intracellular organelles. Endoplasmic reticulum (ER) is set as thin cylinders (radius: 0.08 |mgr|\ m, length 1.8 |mgr|\ m; green colored objects), and a mitochondrion is set as a thich cylinder (radius: 0.2 |mgr|\ m, length 1.8 |mgr|\ m; yellow colored objects). 

|

.. image:: imgs/Scheme.jpg
   :scale: 60%
   :align: center

|

2. The designed shape is embedded in a voxel space for simulation (Lines 1-35, 11_create_dend.py). First, a voxel space of 96 |times| 60 |times| 96 voxels (20 nm/voxel) is set as a numpy 3D array 'vol_dend' (Line 14). Then, the spine head is created using a Python module function morphology.ball (skimage; Line 10; 0.25 |mgr|\ m, 12 voxels), which is added to 'vol_dend' (add_shape; Line 15). Similarly, we make the spine neck and dendrite as cylinders (Lines 11 and 12, respectively), and add them to 'vol_dend' (Lines 16 and 17, respectively). The in-house function 'add_shape,' simply overlays overwrapped regions, and represents the filled areas as 1 (cytosolic region) and the void areas as 0 (extracellular region). Thus, we have already built the contour of the spiny dendrite in 'vol_dend.' Similarly, the PSD, Mitochondrion, and ER are embedded in vol_psd, vol_mito, and vol_er, respectively.

|

3. Here, the volume sizes need to be corrected. LM can only simulate the volume with a size of a multiple of 32 |times| 32 |times| 32 voxels, whereas the current voxel size is 96 |times| 60 |times| 96 voxels. We thus execute a utility function 'lmpad' (Lines 39-42). The voxel sizes are automatically expanded to the multiple of 32 (96 |times| 64 |times| 96 voxels in this case). Users can of course set a multiple of 32 voxels from the beginning.



.. literalinclude:: ../../tutorial/1/11_create_dend.py
   :language: python
   :linenos:
   :caption: 11_create_dend.py


4. Then, we use the CreateSurface class to generate smoothed surfaces of each object (Lines 53, 63, and 67). The generated surfaces are required for locating surface molecules as well as visualizing simulation results. Each surface is composed of triangles that are specified by three vertices (Lines 55 and 56, respectively). Also, the CreateSurface class can generate the surface areas per volume. LD can distribute surface molecules in the voxel space, depending on the surface areas per volume (Line 56). We can further select the surface triangles within the areas of PSD (face_id_psd, Lines 58, 59), to distribute molecules only in this area.

|

5. Finally, the generated variables are assembled in a Python dictionary variable 'm' (Lines 51-77), and 'm' is saved into the HDF container file 'models/ball_and_stick.h5' (Lines 79-85).

.. literalinclude:: ../../tutorial/1/12_show_dend.py
   :language: python
   :linenos:
   :caption: 12_show_dend.py

.. image:: imgs/ball_and_stick.png
   :scale: 50%
   :align: center

6. Execute 'python3 11_create_dend.py'. If users have successfully created the volumes and surfaces, the subsequent script 'python3 12_show_dend.py' will show its 3D shape (Figure above). We use the spiny dendrite for simulation.

