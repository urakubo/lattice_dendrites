def set_molecules(cell):

	print('\nDefine molecules')
	names_concs = {'YFP': 1, 'bleached YFP': 0}
	names = list(names_concs.keys())
	cell.define_species( names )

	print('\nAdd molecules')
	cyt = 'cytosol'
	for name, conc in names_concs.items():
		cell.add_molecule_uM(name, conc, cyt)

	print('\nSet diffusion')
	d_yfp = 0.7 * 1e-12 # (0.7 um2/s) Kang 2012; Traffic 13:1589-1600
	for name in names:
		cell.set_diffusion(name, d_yfp, domain_name=cyt)
