=============
Label a spine
=============

We would like to label a spine volume to obtain molecular concentration of this region. To enable it, HU also developed software, UNI-EM annotator, to manually label any region-of-interests on the surface of objects, and it also serves to obtain volumes within the labeled regions. 

Launch UNI-EM and open the 'annot_ball_and_stick' folder from 'Open Annotator Folder' of the pulldown menu 'Annotator'. Label the spine head in the UNI-EM annotator. Labeled areas are automatically saved.

.. image:: imgs/UNI-EM_annot2.jpg
   :scale: 50%
   :align: center

Convert the labeled areas to labeled volumes, and save them.


.. literalinclude:: ../../tutorial/1/main2_label_spine1.py
   :language: python
   :linenos:


.. image:: imgs/painted.jpg
   :scale: 50%
   :align: center


Confirm the successful segmentation in the voxel space by visualizing it.


.. literalinclude:: ../../tutorial/1/main2_label_spine2.py
   :language: python
   :linenos:


.. image:: imgs/labels_ball_and_stick.png
   :scale: 50%
   :align: center


The spine has a geometric shape. We can thus programmably label the spine volume, not using the UNI-EM annotator, as follows.


.. literalinclude:: ../../tutorial/1/main2_label_spine3.py
   :language: python
   :linenos:


That is all for labeling.
