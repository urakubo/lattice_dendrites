from tut1_functions import show_dendrite

input_morpho_file = "models/ball_and_stick.h5"
output_image_file = "imgs/ball_and_stick.png"

mlab = show_dendrite(input_morpho_file)
mlab.savefig(output_image_file)
mlab.show()