import h5py
from mayavi import mlab

# Load surface meshes
input_morpho_file = "files/ball_and_stick.h5"
output_image_file = "imgs/ball_and_stick.png"

with h5py.File(input_morpho_file,'r') as f:
	bound_v   = f['boundary vertices'][()]
	bound_f   = f['boundary faces'][()]
	PSD_ids   = f['PSD ids in boundary faces'][()]
	mito_v    = f['mitochondrion vertices'][()]
	mito_f    = f['mitochondrion faces'][()]
	er_v      = f['er vertices'][()]
	er_f      = f['er faces'][()]
	pitch     = f['unit length per voxel (um)'][()]

bound_v = bound_v / pitch
mito_v = mito_v / pitch
er_v   = er_v / pitch

# Plot surface mesh
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(700,700))
mlab.view(90, 90, 300, [ 50, 30, 50 ] )

mlab.triangular_mesh(mito_v[:,0] , mito_v[:,1] , mito_v[:,2] , mito_f, color=(1.0,1.0,0.6), opacity=0.6)
mlab.triangular_mesh(er_v[:,0] , er_v[:,1] , er_v[:,2] , er_f, color=(0.7,1.0,0.7), opacity=0.6)
mlab.triangular_mesh(bound_v[:,0], bound_v[:,1], bound_v[:,2], bound_f, color=(0.7,0.7,0.7), opacity=0.3)
mlab.triangular_mesh(bound_v[:,0], bound_v[:,1], bound_v[:,2], bound_f[PSD_ids,:], color=(1,0,0), opacity=0.3)

mlab.savefig(output_image_file)
mlab.show()