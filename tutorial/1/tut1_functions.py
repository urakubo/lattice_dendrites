import numpy as np
import h5py
from skimage import morphology
from mayavi import mlab

def make_cylinder(radius, length, direction = 2):
	arrays = [morphology.disk(radius) for _ in range(length)]
	cylinder = np.stack(arrays, axis=0)
	cylinder = cylinder.swapaxes(0, direction)
	return cylinder

def add_shape(volume, object, loc_center):
	s = np.array(object.shape)
	c = np.floor(s/2).astype(int)
	b = loc_center - c
	e = b + s
	volume[b[0]:e[0], b[1]:e[1], b[2]:e[2] ] += object
	volume = (volume > 0).astype(np.uint8)
	return volume

def show_dendrite(input_morpho_file):

	with h5py.File(input_morpho_file,'r') as f:
		bound_v   = f['bound vertices'][()]
		bound_f   = f['bound faces'][()]
		PSD_ids   = f['face id psd'][()]
		mito_v    = f['mito vertices'][()]
		mito_f    = f['mito faces'][()]
		er_v      = f['er vertices'][()]
		er_f      = f['er faces'][()]
		pitch     = f['unit length (um)'][()]

	bound_v = bound_v / pitch
	mito_v = mito_v / pitch
	er_v   = er_v / pitch

	# Plot surface mesh
	mlab.figure(bgcolor=(1.0,1.0,1.0), size=(700,700))
	mlab.view( 90, 90, 300, [ 50, 30, 50 ] )

	mlab.triangular_mesh(mito_v[:,0] , mito_v[:,1] , mito_v[:,2] , mito_f, color=(1.0,1.0,0.6), opacity=0.6)
	mlab.triangular_mesh(er_v[:,0] , er_v[:,1] , er_v[:,2] , er_f, color=(0.7,1.0,0.7), opacity=0.6)
	mlab.triangular_mesh(bound_v[:,0], bound_v[:,1], bound_v[:,2], bound_f, color=(0.7,0.7,0.7), opacity=0.3)
	mlab.triangular_mesh(bound_v[:,0], bound_v[:,1], bound_v[:,2], bound_f[PSD_ids,:], color=(1,0,0), opacity=0.3)

	return mlab