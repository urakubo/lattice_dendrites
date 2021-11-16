import sys, os
from pyLD import *

annot_folder      = 'annot_ball_and_stick'
object_id         = 1
output_label_file = 'models/labels_ball_and_stick.h5'

c = CreateLabeledVolumeFromUniEM(annot_folder)
c.exec(object_id)
c.save(output_label_file)