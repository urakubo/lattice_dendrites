=============
Label a spine
=============


For analyses, we often need to obtain molecular concentration at a specific region, such as a spine. To enable this, we label the spine volume. It is easy to do it in the case of geometrically shaped regions, because users can label the specific regions by re-defining the same shapes in the same space. In the script (21_label_head.py), we re-define the spine volume alone and saved it in the container 'label volume' of the file 'models/labels_ball_and_stick.h5'. Similarly, the label id (ID: 1) is saved in the container 'label ids.'

.. literalinclude:: ../../tutorial/1/21_label_head.py
   :language: python
   :linenos:
   :caption: 21_label_head.py

If the spine label has been successfully saved, users can confirm it using the script for visualzation '22_show_label.py.' The labeled area appears as a colored part.

.. literalinclude:: ../../tutorial/1/22_show_label.py
   :language: python
   :linenos:
   :caption: 22_show_label.py


.. image:: imgs/labels_ball_and_stick.png
   :scale: 50%
   :align: center


It is another task to label specific regions of a morphologically realistic cell/neuron. We will try it in tutorial 2.
