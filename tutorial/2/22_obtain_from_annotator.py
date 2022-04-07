import sys, os
from pyLD import *

annot_folder      = 'annot_realistic'
object_id         = 1
output_label_file = 'models/labels_realistic.h5'

c = CreateLabelVolumeFromUniEM(annot_folder)
c.create(object_id)
c.save(output_label_file)