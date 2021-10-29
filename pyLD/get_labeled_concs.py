import numpy as np
import h5py
import os, sys
from .utils import get_volume_info

def count_molecules(particles, targ_spine_labels, S):
	num_molecules_spine = {}
	for Targ in S.keys():
		num_molecules_spine[Targ] = [0 for i in range(len(targ_spine_labels))]

	for j in range(particles.shape[3]):
		p = np.ravel(particles[:,:,:,j])
		for i, targ_spine_label in enumerate(targ_spine_labels):
		    pp = p[targ_spine_label]
		    for Targ in S.keys():
		        num_molecules_spine[Targ][i] += np.count_nonzero( pp == S[Targ] ) 
	return num_molecules_spine


def save_labeled_concs(filename, num_molecules_spine, uMs, timepoints, S, ids_spine):
	print('Savefile : ', filename)
	with h5py.File(filename, 'w') as f:
		f.create_dataset('t', data=timepoints)
		f.create_dataset('ids', data=ids_spine)
		f.create_group('number')
		for Targ in S.keys():
			f['number'].create_dataset(Targ, data=num_molecules_spine[Targ])
		f.create_group('conc in uM')
		for Targ in S.keys():
			f['conc in uM'].create_dataset(Targ, data=uMs[Targ])


def get_labeled_concs(lm_filename, labels, output_filename = None, monitor_species = 'Ca'):

	"""Get time series of moleuclar numbers/concentrations within labeled areas from LM simulation result.

	Args:
		lm_filename (str): Filename of lm simulation.
		labels (numpy[int/bool]): Label volume (3D array). id=0 will be ignored.
		output_filename (none/str): If None, get_labeled_concs will give returns below. If specified, it will save a hdf5 file.
		monitor_species (none/str): If None, get_labeled_concs will give no messages in the console. If a moleuclar species is specified, it will show a example result.

	Returns:
		(tuple): Tuple containing:

		- num_molecules (dict): Time series of the numbers of molecules of the specified molecular species. The dict container has {'species1': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], 'Species2': [[num_label1_t1, num_label2_t1, ...], [num_label1_t2, num_label2_t2, ...], ...], ... }.
		- uMs (dict): Time series of molecular concentrations. The dict container has {'species1': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], 'Species2': [[conc_label1_t1, conc_label2_t1, ...], [conc_label1_t2, conc_label2_t2, ...], ...], ... }.
		- timepoints (numpy[float]): Timepoints (s)
		- ids_spine (numpy[int]): Labels in the volume. The numpy array has [label1, label2, ...].

	"""

    # Decode info
	if not isinstance(lm_filename, str) :
		print("input lm_filename must be str.")
		return False

	cyt = 1
	NA  = 6.022e23
	num_voxels, volume_in_L, spacing, S = get_volume_info(lm_filename, cyt)

	# Decode labels
	ids_spine, nums_spine_voxels = np.unique(labels, return_counts=True)
	ids_spine          = ids_spine[1:] #  0 was removed.
	nums_spine_voxels  = nums_spine_voxels[1:]
	vols_spine_in_L    = nums_spine_voxels * spacing * spacing * spacing * 1000

	targ_spine_labels = []
	labels_flat = np.ravel(labels)
	for id_spine in ids_spine:
	    targ_spine_label = np.where( (labels_flat == id_spine) )
	    targ_spine_labels.append( targ_spine_label )

	# Time frames
	with h5py.File(lm_filename,'r') as file:
		timepoints = file['Simulations']['0000001']['LatticeTimes'][()]
		frames = [key for key in file['Simulations']['0000001']['Lattice'].keys()]
	frames.sort()

	# Monitor messages
	if monitor_species != None:
	    print('file :', input_lm_file)
	    print('Label ids : ', ids_spine)
	    print('Species    : ', ', '.join(list(S.keys())))


	# Get molecules in labeled volumes at each timepoint
	num_molecules = {}
	for Targ in S.keys():
	    num_molecules[Targ] = []

	for f in frames:
		#
		num_molecules_time_i = {}
		for Targ in S.keys():
		    num_molecules_time_i[Targ] = [0 for i in range(max(ids_spine)+1)]

		#
		with h5py.File(lm_filename,'r') as file:
		    particles = file['Simulations']['0000001']['Lattice'][f][:,:,:,:]

		#
		num_molecules_time_i = _count_molecules(particles, targ_spine_labels, S)
		if isinstance(monitor_species, str) and monitor_species in S.keys():
		    print('Monitor species:', monitor_species,', Number_at_time_i:', num_molecules_time_i[monitor_species] )

		#
		for Targ in S.keys():
			num_molecules[Targ].append(num_molecules_time_i[Targ])

	uMs = {}
	for Targ in S.keys():
		num_molecules[Targ] = np.array(num_molecules[Targ])
		uMs[Targ] = num_molecules[Targ] / NA * 1e6 / vols_spine_in_L

	if output_filename == None:
		return num_molecules, uMs, timepoints, ids_spine
	elif isinstance(output_filename, str):
		return _save_labeled_concs(output_filename, num_molecules, uMs, timepoints, S, ids_spine)
	else :
		print('output_filename is not str.')
		return False

if __name__ == '__main__':
    
	input_lm_file    = 'lms900_2_4/_stimrun_00.lm'
	#
	input_label_file = "lm_annot/labels.hdf5"
	with h5py.File(input_label_file,'r') as f:
	    labels = f['dendrite'][()]
	output_file = 'test.hdf5'

	get_labeled_concs(input_lm_file, labels,  output_filename = output_file, monitor_species = 'Ca')
