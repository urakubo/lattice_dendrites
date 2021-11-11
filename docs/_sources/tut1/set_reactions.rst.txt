=============
Set reactions
=============

まずは、好きなディレクトリに移動して、sphinx-quickstartを実行します。

以下では report というドキュメント名を指定して実行しています。


.. code-block:: python
	:linenos:

	import h5py
	import numpy as np

	from pyLM import *
	from pyLM.units import *
	from pyLD import * 

	from setMolecules import setMolecules

	input_filename_morph='morph/ball_and_stick.h5'
	output_filename_lm='morph/ball_and_stick.lm'

	ext  = 'default'
	cyt  = 'cytoplasm'
	er   = 'er'
	mito = 'mito'
	domains = {ext: 0, cyt: 1, er: 2, mito: 3}


	print('\nLoad geometry data.\n')

	with h5py.File(input_filename_morph,'r') as r:
	    dendrite							= r['dendrite'][()]
	    mitochondrion						= r['mitochondrion'][()]
	    dendrite_not_mitochondrion_not_ER	= r['dendrite not mitochondrion not ER'][()]
	    PSD									= r['PSD'][()]
	    vol_bound							= r['boundary areas in volume'][()]
	    pitch								= r['unit length per voxel (um)'][()]
	    ER_area								= r['ER areas in volume'][()]

	volume = (dendrite_not_mitochondrion_not_ER > 0) * domains[cyt]
	volume[mitochondrion > 0] = domains[mito]

	vol_PSD  = vol_bound * (PSD > 0)
	surfaces = {'PSD': vol_PSD, 'cell boundary': vol_bound}


	print('\nSimulation setup.\n')

	latticeSpacing=micron(pitch)
	nx,ny,nz = volume.shape
	sim  = RDME.RDMESimulation(dimensions=micron(nx*pitch,ny*pitch,nz*pitch), spacing=latticeSpacing)
	cell = BuildAnyShape(sim, volume, domains, surfaces)


	NA = 6.022e23
	number_1umol = NA /(1e6)

	num_voxel_cyt = np.count_nonzero(volume == domains[cyt])
	cyt_vol_in_L  = num_voxel_cyt * latticeSpacing * latticeSpacing * latticeSpacing * 1000
	number_1uM    = number_1umol * cyt_vol_in_L

	print('Cytosolic volume in fL:', cyt_vol_in_L * 1e15)
	print('Number of molecules per 1 uM:', number_1uM)

	print('\nSet molecules.\n')

	molecules = setMolecules(sim, cyt)
	molecules.setDiffusion()
	molecules.setReactions()

	conc_Ca  = 100 # uM
	conc_CaM = 100 # uM
	conc_CB  = 120 # uM
	conc_CN  = 5   # uM
	num_Ca   = int(conc_Ca  * number_1uM /2)
	num_CaM  = int(conc_CaM * number_1uM /2) 
	num_CB   = int(conc_CB  * number_1uM /2)
	num_CN   = int(conc_CN  * number_1uM /2)

	print('num_Ca : ', num_Ca )
	print('num_CaM: ', num_CaM)
	print('num_CB : ', num_CB )
	print('num_CN : ', num_CN )

	cell.add_cytosolic_molecules('Ca'  , num_Ca  , cyt) # Absolute number
	cell.add_cytosolic_molecules('N0C0', num_CaM , cyt)
	cell.add_cytosolic_molecules('CB'  , num_CB  , cyt)
	cell.add_cytosolic_molecules('CN'  , num_CN  , cyt)
	# 6000 moleules per 100uM CaM and 0.1fL Spine

	# 488e-12 [mol/m2] => (6.02e23) * 488e-12 / 1e12  [number/um2] 

	PMCA  =  488e-12*NA/1e12 # number per um2
	NCX   =  488e-12*NA/1e12 # number per um2
	NMDAR =  500e-12*NA/1e12 # number per um2
	RyR   =  488e-12*NA/1e12 # number per um2

	cell.add_surface_molecules('PMCA', PMCA/2 , 'cell boundary')
	cell.add_surface_molecules('NCX' , NCX/2  , 'cell boundary')
	cell.add_surface_molecules('NR'  , NMDAR/2, 'PSD')
	#cell.add_er_molecules('RC1c0' , RyR)

	print('\nSet time.\n')
	sim.setTimestep(microsecond(3.0))
	sim.setWriteInterval(0.05)
	sim.setLatticeWriteInterval(0.05)
	sim.setSimulationTime(5.0)

	print('\nSave simulation setup.\n')
	if(os.path.isfile(output_filename_lm) == True):
	    os.remove(output_filename_lm)
	sim.save(output_filename_lm)

	for molecular_name in molecules.cyt_mol:
	    particleNum=sim.particleMap[molecular_name]
	    print(molecular_name, str(particleNum))
	    
	for molecular_name in molecules.sur_mol:
	    particleNum=sim.particleMap[molecular_name]
	    print(molecular_name, str(particleNum))



That is all for simulation setup.
