=============
Label a spine
=============

We would like to label a spine volume to obtain molecular concentration of this region.
Because the spine has a geometric shape, we can programmably label the spine volume, not using the UNI-EM annotator, as follows.

Convert the labeled areas to labeled volumes, and save them.


.. literalinclude:: ../../tutorial/1/21_label_head.py
   :language: python
   :linenos:
   :caption: 21_label_head.py


Confirm the successful segmentation in the voxel space by visualizing it.


.. literalinclude:: ../../tutorial/1/22_show_label.py
   :language: python
   :linenos:
   :caption: 22_show_label.py


.. image:: imgs/labels_ball_and_stick.png
   :scale: 50%
   :align: center




That is all for labeling.
