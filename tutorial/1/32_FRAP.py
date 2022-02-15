def set_molecules(cell):

	print('\nDefine molecules\n')
	name_yfp = 'YFP'
	cell.define_species(name_yfp)

	print('\nAdd molecules\n')
	conc_yfp = 1 # uM
	cyt      = 'cytosol'
	cell.add_molecule_uM(name_yfp, conc_yfp, domain_name=cyt)

	print('\nSet diffusion\n')
	d_yfp    = 0.7 * 1e-12 # (0.7 um2/s) Kang 2012; Traffic 13:1589-1600
	cell.set_diffusion(name_yfp, d_yfp, domain_name=cyt)