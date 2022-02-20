import os
from tut2_functions import show_dendrite

input_morpho_file = "models/realistic_dendrite.h5"
output_image_file = "imgs/ball_and_stick.png"
os.makedirs("imgs", exist_ok=True)

mlab = show_dendrite(input_morpho_file)
mlab.savefig(output_image_file)
#mlab.show()