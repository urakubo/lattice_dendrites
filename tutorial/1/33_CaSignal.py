def set_molecules(cell):

	print('\nDefine molecules\n')
	name_Ca        = 'Ca'
	surf_numbers   = {'NMDAR_act':0    ,'NMDAR_inact':20   , 'Pump':10000, 'Pump_Ca':0}
	surf_local     = {'NMDAR_act':'PSD','NMDAR_inact':'PSD', \
			'Pump':'cell boundary', 'CaPump':'cell boundary'}
	cell.define_species( name_Ca )
	cell.define_species( list(surf_numbers.keys()) )

	print('\nAdd molecules\n')
	conc_Ca        = 0
	cyt            = 'cytosol'
	cell.add_molecule(name_Ca, conc_Ca, cyt)
	for k, v in surf_numbers.items():
		cell.add_surface_molecule(k, v, surf_local[k])

	print('\nSet diffusion\n')
	diff_Ca        = 10 * 1e-12
	cell.set_diffusion(name_Ca, diff_Ca, cyt)

	print('\nSet reactions\n')
	k_nmdar_deact = 0
	k_channel_Ca  = 0
	k_pump_Ca     = 0
	kon_Ca        = 0
	koff_Ca       = 0
	oneway_reacs = \
		[['NMDAR_act',       'NMDAR_inact' , k_nmdar_deact], \
		 ['NMDAR_act', ('NMDAR_act', 'Ca') , k_channel_Ca] , \
		 ['Pump_Ca'  ,              'Pump' , k_pump_Ca]]
	for p in oneway_reacs:
		cell.reac_oneway_uM(reac=p[0], prod=p[1], rate=p[2], domain_name=cyt)
	#
	cell.reac_twoway_uM(reac=('Pump', 'Ca'), \
		prod='Pump_Ca', \
		rates=(kon_Ca, koff_Ca), \
		domain_name=cyt)