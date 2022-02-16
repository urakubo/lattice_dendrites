def set_molecules(cell):

	print('\nDefine molecules')
	name_Ca        = 'Ca'
	surf_numbers   = {'active NMDAR':0    ,'inactive NMDAR':300, \
			'Pump':100, 'Pump-Ca':0}
	surf_local     = {'active NMDAR':'PSD','inactive NMDAR':'PSD', \
			'Pump':'cell boundary', 'Pump-Ca':'cell boundary'}
	cell.define_species( name_Ca )
	cell.define_species( list(surf_numbers.keys()) )

	print('\nAdd molecules')
	conc_Ca        = 0
	cyt            = 'cytosol'
	cell.add_molecule(name_Ca, conc_Ca, cyt)
	for k, v in surf_numbers.items():
		cell.add_surface_molecule(k, v, surf_local[k])

	print('\nSet diffusion')
	diff_Ca        = 10 * 1e-12
	cell.set_diffusion(name_Ca, diff_Ca, cyt)

	print('\nSet reactions')
	k_nmdar_deact = 20
	k_channel_Ca  = 1000
	k_pump_Ca     = 10
	kon_Ca        = 10
	koff_Ca       = 1
	oneway_reacs = \
		[['active NMDAR', 'inactive NMDAR' , k_nmdar_deact], \
		 ['active NMDAR', ('active NMDAR', 'Ca') , k_channel_Ca] , \
		 ['Pump-Ca'     , 'Pump' , k_pump_Ca]]
	for p in oneway_reacs:
		cell.reac_oneway_uM(reac=p[0], prod=p[1], rate=p[2], domain_name=cyt)
	#
	cell.reac_twoway_uM(reac=('Pump', 'Ca'), \
		prod='Pump-Ca', \
		rates=(kon_Ca, koff_Ca), \
		domain_name=cyt)
