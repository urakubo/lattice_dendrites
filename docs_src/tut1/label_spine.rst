=============
Label a spine
=============


For analyses, we often need to obtain molecular concentration of a specific region, such as a spine. To enable this, we label a spine volume. It is easy to do it in the case of the geometrically shaped cells, because users can label a specific area by re-defining the same sphere in the same space. In the sample script (12_label_head.py), the total cytosolic volume 'vol_dend_not_mito_not_er' is loaded from 'models/ball_and_stick.h5', and saved as the container 'ref volume' in the label file 'models/labels_ball_and_stick.h5'. Similarly, the label id is saved as the container 'label ids'.


.. literalinclude:: ../../tutorial/1/21_label_head.py
   :language: python
   :linenos:
   :caption: 21_label_head.py

If users can successfully save the spine label, they can confirm it using the visualzation script ‘22_show_label.py’. The labeled area appears as a colored part.

.. literalinclude:: ../../tutorial/1/22_show_label.py
   :language: python
   :linenos:
   :caption: 22_show_label.py


.. image:: imgs/labels_ball_and_stick.png
   :scale: 50%
   :align: center


It is not easy to label specific regions in the case of a morphologically realistic model. We will try it in tutorial 2. 
