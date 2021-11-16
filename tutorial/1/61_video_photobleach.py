import numpy as np
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


def remove(lattice, sys_param, event_param):
    i    = sys_param['i']
    time = sys_param['time']
    print('\nRemove event at: {:g}, Current time: {:.3f}\n'.format(i, time))

    target_volume = (sys_param['label volume'] == event_param['target label id'] )
    tmp_lattice = lattice
    for i in range(lattice.shape[3]):
        tmp_lattice = lattice[:,:,:,i]
        tmp_lattice[target_volume] = 0
        lattice[:,:,:,i] = tmp_lattice

    NR   = []
    for i in range(particles.shape[3]):
        NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR']).tolist() )
        NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR_O']).tolist() )
        NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR_Glu']).tolist() )
        NR = np.unravel_index(NR, particles[:,:,:,0].shape )
        plot_points2 = mlab.points3d(NR[0], NR[1], NR[2], color=(1,0,0), scale_factor=1.5,line_width=0.1)

        for i in range(particles.shape[3]):
            for id in target_molecules:
                Molecule.extend( np.flatnonzero(particles[:,:,:,i] == S[id] ).tolist() )
        if  Molecule != []: 
            M  = np.unravel_index(Molecule, particles[:,:,:,0].shape )
            plot_points1 = mlab.points3d(M[0], M[1], M[2], color=col, scale_factor=4.0,line_width=0.1)
  
        text = mlab.text(0.1,0.05,"{0:.3f} s".format(t) ,color=(0,0,0), width=0.3)


    return event_param



print('mlab.options.offscreen: ', mlab.options.offscreen)
image_id    = 0
for lmfile in input_lm_files:
    print('file :', lmfile)
    hfile = h5py.File(lmfile,'r')
    Timepoints = hfile['Simulations']['0000001']['LatticeTimes'][()]
    Timepoints = Timepoints + Time_offset
    Timepoints = Timepoints.tolist()
    Timepoints = Timepoints[1:]
    Time_offset = Timepoints[-1]
    
    Frames = list(hfile['Simulations']['0000001']['Lattice'].keys())
    Frames.sort()
    Frames.pop(0)
    
    for t, f in zip(Timepoints, Frames):
        print('Timepoint: ', t)
        if (t < -2):
            continue
        if (t > 5):
            continue
        particles = hfile['Simulations']['0000001']['Lattice'][f][:,:,:,:]
        Molecule = []
        NR   = []
        for i in range(particles.shape[3]):
            NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR']).tolist() )
            NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR_O']).tolist() )
            NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR_Glu']).tolist() )
        NR = np.unravel_index(NR, particles[:,:,:,0].shape )
        plot_points2 = mlab.points3d(NR[0], NR[1], NR[2], color=(1,0,0), scale_factor=1.5,line_width=0.1)
        for i in range(particles.shape[3]):
            for id in target_molecules:
                Molecule.extend( np.flatnonzero(particles[:,:,:,i] == S[id] ).tolist() )
        if  Molecule != []: 
            M  = np.unravel_index(Molecule, particles[:,:,:,0].shape )
            plot_points1 = mlab.points3d(M[0], M[1], M[2], color=col, scale_factor=4.0,line_width=0.1)
  
        text = mlab.text(0.1,0.05,"{0:.3f} s".format(t) ,color=(0,0,0), width=0.3)
        # #print(dir(text.property.font_size))
        #text.property.font_size = 10
        #text.property.font_family = 'arial'
        ##mlab.view(-180.0-image_id*0.5, 90, 2000, [120.0, 120.0, 480.0 ])
        mlab.view(-180.0, 90, 2000, [120.0, 120.0, 480.0 ])
        ffname = join(output_dir, output_file +str(image_id).zfill(4)+ '.png')
        mlab.savefig(ffname)
        if  Molecule != []: 
            plot_points1.remove()
        plot_points2.remove()
        text.remove()
        image_id = image_id + 1
    hfile.close()

##
## Generate video
##
ffname = join(output_dir, output_file +'%04d'+ '.png')
com = ['ffmpeg','-r', '10', '-i', ffname,'-pix_fmt', 'yuv420p', target_molecules[0]+suffix]
print(' '.join(com))
s.call(com)

# ffmpeg -r 10 -i png_Ca/test_%04d.png -pix_fmt yuv420p Ca_4.mp

