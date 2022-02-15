import os
from tut1_functions import show_dendrite

input_morpho_file = "models/ball_and_stick.h5"
output_image_file = "imgs/ball_and_stick.png"
os.makedirs("imgs", exist_ok=True)

mlab = show_dendrite(input_morpho_file)
mlab.savefig(output_image_file)
mlab.show()