from __future__ import print_function
from __future__ import division

import sys, os, time, errno
from os import path, pardir
import numpy as np
import h5py
import json
import pathlib
import trimesh



def save_uniem_annotator(foldername, pitch, ids_volume, surfaces):

	"""Save a volume and surface meshes under the style of UNI-EM annotator.
	This is aimed to directly connect pyLD-generated surfaces with UNI-EM annotator.

	Args:
		foldername (str): Relative/absolute path to a foldername for uni-em annotator. The folder will be generated if it does not exit.
	    pitch (float): xyz pitch per voxel (um)
	    volume (numpy[int/bool]): Input volume (3D array)
	    surfaces (dict): Surface meshes. The dict container must have {id1: [vertices1, faces1, col1], id2: [vertices2, faces2, col2], ... } where id: int,  vertices: 3xX float np.array,  faces: 3xY int np.array, and col: [R(0-255), G(0-255), B(0-255)].

	Returns:
		(bool): True or false
	"""
	# Initialization
	if not isinstance(surfaces, dict):
		print('surfaces must be dict, but ', type(surfaces))
		return False

	if ids_volume.ndim != 3 :
		print('ids_volume must be numpy 3D array, but ', ids_volume.ndime)
		return False

	if ids_volume.dtype not in ['int8','int16','int32','int64','uint8','uint16','uint64']:
		print('ids_volume must be numpy integers, but ', ids_volume.dtype)
		return False

	if not isinstance(foldername, str):
		print('foldername must be str, but ', type(foldername))
		return False

	p = pathlib.Path(foldername)
	if p.is_absolute() == False:
		folder = p.resolve()
	else:
		folder = foldername

	params = Params(folder)

	# Generation
	if os.path.isdir(folder) == False:
		os.makedirs(folder)
	save_uniem_annotator_preprocess(params)
	save_uniem_annotator_generate_info_file(ids_volume, params.surfaces_segment_info_json_file, surfaces)
	save_uniem_annotator_mesh(params, surfaces)
	save_uniem_annotator_post_process(params, pitch, ids_volume)
	return True 


def save_uniem_annotator_mesh(params, surfaces):
	for id, vf in surfaces.items():
		mesh = trimesh.Trimesh(vertices=vf[0], faces=vf[1])
		filename = os.path.join(params.surfaces_whole_path ,str(id).zfill(10)+'.stl')
		mesh.export(file_obj=filename)
	return True


def save_uniem_annotator_preprocess(params):

	print('Annotator folder is being generated.')
	targ = params

	if os.path.isdir(targ.surfaces_path) == False:
		os.makedirs(targ.surfaces_path)
		os.makedirs(targ.surfaces_whole_path)
	if os.path.isdir(targ.skeletons_path) == False:
		os.makedirs(targ.skeletons_path)
		os.makedirs(targ.skeletons_whole_path)
	if os.path.isdir(targ.volume_path) == False:
		os.makedirs(targ.volume_path)
	if os.path.isdir(targ.paint_path) == False:
		os.makedirs(targ.paint_path)

	return True


def save_uniem_annotator_post_process(params, pitch, ids_volume):

    targ = params
	##
    ph = pitch
    pw = pitch
    pz = pitch

    wmax = ids_volume.shape[0]
    hmax = ids_volume.shape[1]
    zmax = ids_volume.shape[2]
	##
    with h5py.File(targ.volume_file, 'w') as f:
    	f.create_dataset('volume', data=ids_volume)
	##
    data_dict = {
		'boundingbox_voxel':{
			'x': wmax,
			'y': hmax,
			'z': zmax
			},
		'boundingbox_um':{
			'x': pw * wmax,
			'y': ph * hmax,
			'z': pz * zmax
			},
		'pitch_um':{
			'x': pw,
			'y': ph,
			'z': pz
			},
		}
    with open( targ.surfaces_volume_description_json_file , 'w') as f:
    	json.dump(data_dict, f, indent=2, ensure_ascii=False)

    print('')
    print('Annotator folder generated.')
    print('')

    return True


def save_uniem_annotator_generate_info_file(ids_volume, surfaces_segment_info_json_file, surfaces):
    ids_nums = np.unique(ids_volume, return_counts=True)
    ids   = ids_nums[0]
    names = [str(id).zfill(10) for id in ids]
    sizes = ids_nums[1]
    colormap = np.random.randint(255, size=(ids.shape[0], 3), dtype='int')

    for id, vf in surfaces.items():
    	colormap[int(id), 0] = vf[2][0]
    	colormap[int(id), 1] = vf[2][1]
    	colormap[int(id), 2] = vf[2][2]

    if ids[0] == 0:
    	ids   = np.delete(ids, 0)
    	names.pop(0)
    	sizes = np.delete(sizes, 0)
    	colormap = np.delete(colormap, 0, 0)

    ids      = ids.tolist()
    sizes    = sizes.tolist()
    colormap = colormap.tolist()

    print('Constainer shape: ', ids_volume.shape)
    print('IDs  : ', ids)
    print('names: ', names)
    print('sizes: ', sizes)
    print('cols : ', colormap)

	##
    keys = ['id', 'name', 'size']
    data_dict = [dict(zip(keys, valuerecord)) for valuerecord in zip(ids, names, sizes)]

    for i in range(len(data_dict)):
    	col = {'confidence': 0, 'r': colormap[i][0], 'g': colormap[i][1],  'b': colormap[i][2],  'act': 0}
    	data_dict[i].update(col)

    print('data_dict: ', data_dict)

    with open( surfaces_segment_info_json_file , 'w') as f:
    	json.dump(data_dict, f, indent=2, ensure_ascii=False)

    return True



