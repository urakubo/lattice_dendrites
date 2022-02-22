.. include:: ../isonum.txt
.. include:: ../isogrk1.txt

=========================
Label spines using UNI-EM
=========================

In tutorial 1, the target spine was labeled by its re-generation. This method cannot be used for morphologically realistic spines. We thus developed the software to manually label the realistic spines or any other region-of-interests, which is named as UNI-EM annotator. 

The installation and usage of UNI-EM annotator are described elsewhere. We here introduce the function save_uniem_annotator to convert data to the format of UNI-EM annotator (Line 26 in 21_convert_to_annotator.py). In this function, 'pitch' denotes the unit length of each voxel, and 'volume' contains the objects of the realistic dendrite (1: cytosol, 2: mitochondrion, 3: ER).The dict variable 'surfaces' contains the volume ids, surface vertices and faces, and colors. This function generates the files of UNI-EM annotator in the directory specified by 'annot_folder'.


.. literalinclude:: ../../tutorial/2/21_convert_to_annotator.py
   :language: python
   :linenos:
   :caption: 21_convert_to_annotator.py

|

UNI-EM annotator has a paint function to label 3D surfaces. Users can label any region-of-interests (Figure below).


.. image:: imgs/UNI-EM.png
   :scale: 50%
   :align: center

|

The class method 'exec' generates the volume that contains labeled volumes. This function is realized by the Python module Pymeshfix that makes closed surface meshes and re-voxelizes the closed surfaces (Line 9; Figure below). The label volume is saved in the container 'label volume' of the HDF file 'models/labels_realistic.h5' (Lines 6, 10).


.. literalinclude:: ../../tutorial/2/22_obtain_from_annotator.py
   :language: python
   :linenos:
   :caption: 22_obtain_from_annotator.py


.. image:: imgs/spines.png
   :scale: 50%
   :align: center

|

Finally, the labeled surface regions are visualized by the script '23_show_label.py'. The labeled spines would be colored as shown in Figure below. Based on the labeled shape, users can simulate molecular interactions and evaluate the simulation results.

.. literalinclude:: ../../tutorial/2/23_show_label.py
   :language: python
   :linenos:
   :caption: 23_show_label.py


.. image:: imgs/labels_realistic.png
   :scale: 100%
   :align: center


