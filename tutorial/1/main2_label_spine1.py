
import sys, os
main_dir = os.path.abspath(os.path.dirname(sys.argv[0]))  # Dir of main
sys.path.append(os.path.join(main_dir, '..', 'github_pages'))
from pyLD import *

annot_folder    = 'annot_ball_and_stick'
object_id = 1
output_label_file = 'labels_ball_and_stick.h5'

c = CreateLabeledVolumeFromUniEM(annot_folder)
c.exec(object_id)
c.save(output_label_file)



