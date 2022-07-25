# from numbalsoda import lsoda_sig, lsoda
# from numba import njit, cfunc
# import numba as nb
from numba import njit, jit, float64, int64, f8, i8, b1, void, literal_unroll
from numba.typed import Dict, List
from numba.types import types


import os, sys, glob
import numpy as np
import itertools
from scipy.integrate import solve_ivp
import networkx as nx
from MembraneMechanisms import DistributedNa, DistributedK, PointCurrent, PointNMDAR

# from pathlib import Path
import h5py
import pprint
import pickle
import trimesh




# @staticmethod
@njit(nopython=True)
def _ff(t, u,
	cm,
	rm,
	ri_inv,
	connector_to_compartments0,
	connector_to_compartments1,
	compartment_to_connectors0,
	compartment_to_connectors1,
	n_compartments,
	v_connectors,
	ud_num
	):
	
	El = -0.070
	
	# Variables:
	# v : membrane potentials at compartments
	# ud: internal states (distributed mechanisms)
	# up: internal states (point mechanisms)
	# print('t ', t)
	v, ud, up = np.split(u, [n_compartments, ud_num])
	ud        = ud.reshape((-1, n_compartments))
	dv        = np.zeros_like(v)
	# Passive channel
	dv += - (v-El) / rm # passive property

	# Voltage transfer among compartments
	#v_connectors = np.zeros(n_connectors, dtype=float64)
	#for id_con, id_comps in connector_to_compartments.items():
	for id_con, id_comps in zip(connector_to_compartments0, connector_to_compartments1) :
		v_connectors[id_con] = 0
		denominator = 0
		numerator   = 0
		for id_comp in id_comps:
			denominator += ri_inv[id_comp]
			numerator   += v[id_comp] * ri_inv[id_comp]
		v_connectors[id_con] = numerator / denominator
	for id_comp, id_cons in zip(compartment_to_connectors0, compartment_to_connectors1):
		current = 0
		for id_con in id_cons:
			current  += ri_inv[id_comp] * ( v_connectors[id_con] - v[id_comp] )
		dv[id_comp]  += current

	'''
	# Distributed mechanisms
	i_s = 0
	d_ud = np.empty_like(ud)
	for mechanism, num_variables in zip(self.distributed_mechanisms, self.distributed_num_variables):
		
		Im, dd_ud = mechanism(v, dv, ud[i_s:i_s+num_variables, :])
		d_ud[i_s:i_s+num_variables, :] = dd_ud
		dv  += Im
		i_s += num_variables

	# Point mechanisms
	i_s = 0
	I = np.zeros_like(v)
	d_up = np.empty_like(up)
	for mechanism, param, num_variables, loc in \
		zip(self.point_mechanisms, self.point_params, self.point_num_variables, self.point_locations):
		I, d_up[i_s:i_s+num_variables] = mechanism(t, v[loc], up[i_s:i_s+num_variables], param)
		i_s += num_variables
		dv[loc] += I
	'''

	# Altogether
	dv   = dv / cm
	#dudt = np.vstack((dv, d_ud)).reshape(-1)
	#dudt = np.hstack([dv, d_up])

	return dv



class SimulateMembranePotential(DistributedNa, DistributedK, PointNMDAR, PointCurrent):
	"""Simulate membrane potential of a multicompartment model. This was purely written in Python to simulate a multicompartment model that has a realistic shape. Users can insert distributed and point mechanisms. The distributed mechanisms denotes the distributed channels of Na, K, and Ca. Only one instance for each mechanism can be generated. The point mechanisms denotes the currents via NMDARs and patch pipettes. Multple point mechanisms can be inserted simultaneously.

	Args:
		(None):

	Returns:
		(pyLD.SimulateMembranePotential): SimulateMembranePotential object that has the follwing instances:
		- init_V (float): Initial membrane potential
		- p (dict): Parameters of pasive properties 'Cm' (uF/um2), 'Rm' (MOhm-um2), 'El' (V), and 'Ra' (MOhm-um).
		- n_compartments (int): Number of compartments.
		- s (numpy[float]): Surface areas of compartments (1D array; [area1, area2, area3, ...]).
		- ra_inv (numpy[float]): Inverse of axial resistances between compartments (1D array; [1?R1, 1/R2, 1/R3, ...]).
		- tstart (float): Start time of simulation.
		- tend (float): End time of simulation.
	"""

	def __init__(self):
		# Passive properties
		self.init_V         = -0.07
		self.p              = {'Cm':1.0e-8, 'Rm': 250e3, 'El': self.init_V, 'Ra':1.5}
		# Units: Cm = uF/um2,  Rm = MOhm-um2, El = Volt, Ra = MOhm-um
		# Other units: Time = sec, Membrane potentials = Volt


		# Distributed mechanisms
		self.distributed_mechanisms    = []
		self.distributed_num_variables = []
		self.distributed_init          = []
		# Point mechanisms
		self.point_mechanisms          = []
		self.point_num_variables       = []
		self.point_init_variables      = []
		self.point_init                = []
		self.point_event               = []
		self.point_locations           = []
		self.point_params              = []

		self.tstart = 0.0
		self.tend   = 50.0

		PointNMDAR.__init__(self)
		PointCurrent.__init__(self)



	def load_spatial_model_graph(self, graph ):
		"""Import the model described using a networkx graph. Here, each edge correpond to a compartment, and each node is their connector.
		Args:
			graph (networkx.classes.graph): Connections to connectors
		"""
		edges_for_remove = []
		if type(graph).__module__ != 'networkx.classes.graph':
			raise ValueError( 'graph is not "networkx.classes.graph".' )
		else :
			for id, edge in graph.edges.items():
				#print('edge: ', edge)
				#print(" ('length' in edge) ", ('length' in edge))
				if  not 'length' in edge:
					raise ValueError( 'edge {} does not have "length".'.format(id) )
				elif not 'area' in edge:
					raise ValueError( 'edge {} does not have "area".'.format(id) )
				elif not 'volume' in edge:
					raise ValueError( 'edge {} does not have "volume".'.format(id) )
				elif edge['length'] <= 0:
					print('This edge has zero-or-minus length: ', edge)
					print('Removed.')
					edges_for_remove.append(id)
				elif edge['volume'] <= 0:
					print('This edge has zero-or-minus volume: ', edge)
					print('Removed.')
					edges_for_remove.append(id)

		self.graph = graph # Networkx where the edges contain "surface_area", "volume", "length"
		self.graph.remove_edges_from(edges_for_remove)
		
		
		# 
		self.n_compartments = self.graph.number_of_edges()
		print('self.n_compartments ',  self.n_compartments)
		
		self.rm     = np.zeros(self.n_compartments, dtype=float)
		self.ri_inv = np.zeros(self.n_compartments, dtype=float)
		self.cm     = np.zeros(self.n_compartments, dtype=float)
		
		id_model = 0
		for edge in self.graph.edges.values():
			length = edge['length']
			area   = edge['area']
			volume = edge['volume']
			section_area = volume / length
			ri                    = self.p['Ra'] * length / section_area
			self.ri_inv[id_model] = 0.5 / ri
			self.rm[id_model]     = self.p['Rm'] * area
			self.cm[id_model]     = self.p['Cm'] * area
			edge['id_model']      = id_model
			id_model += 1
		
		#print('self.ri_inv ', self.ri_inv)
		#print('self.rm     ', self.rm)
		#print('self.cm     ', self.cm)
		
		
		self.ids_compartment =  [ e['id_model'] for e in graph.edges.values()]
		self.ids_connectors  = sorted( set(itertools.chain(*graph.edges.keys())) )
		# print( 'self.ids_connectors  ', self.ids_connectors )
		self.n_connectors    = max(self.ids_connectors)+1

		self.compartment_to_connectors = {e['id_model']: list(con) for con, e in graph.edges.items()}
		tmp = {node: list(graph.edges(node, data='id_model')) for node in graph.nodes.keys()}
		self.connector_to_compartments = { node: [v[2] for v in values] for node, values in tmp.items() }
		self.connector_to_compartments = { node: values for node, values in self.connector_to_compartments.items() if values != [] }



	def solve_ivp_numba(self, tspan, u0, method='LSODA'):

		compartment_to_connectors0 = List()
		compartment_to_connectors1 = List()
		for  key, value in self.compartment_to_connectors.items():
			compartment_to_connectors0.append(key)
			compartment_to_connectors1.append(value)

		connector_to_compartments0 = List()
		connector_to_compartments1 = List()
		for key, value in self.connector_to_compartments.items():
			connector_to_compartments0.append(key)
			connector_to_compartments1.append(value)

		connector_to_compartments1

		v_connectors = np.zeros(self.n_connectors, dtype=float)


		print('u0.shape ' ,u0.shape)
		print('n_compartments ' ,self.n_compartments)
		print('ud_num ' ,self.ud_num)


		sol = solve_ivp(_ff, tspan, u0, method=method, \
			args=(	self.cm,
					self.rm,
					self.ri_inv,
					connector_to_compartments0,
					connector_to_compartments1,
					compartment_to_connectors0,
					compartment_to_connectors1,
					self.n_compartments,
					v_connectors,
					self.ud_num) )

		return sol

	def run(self):
		"""Run the simulation.

		Args:
			None

		Returns: bool 
			(bool): True if succeeded. Also simulation results are stored in the follwing instances:
			- t (numpy[float]): Time (1D array; [t1, t2, t3, ...]).
			- y (numpy[float]): Variables (2D array; [[v1(t=t1), v2(t=t1), v3(t=t1), ...], [v1(t=t2), v2(t=t2), v3(t=t2), ...], [v1(t=t3), v2(t=t3), v3(t=t3), ...], ...]
		"""

		# Init distributed channels
		ud0 = [self.init_V]
		for init in self.distributed_init:
			ud0.extend( init() )
		ud0 = np.tile(np.array(ud0),(self.n_compartments,1)).T.reshape(-1)
		ud_num = ud0.shape[0]
		self.ud_num = ud_num

		# Init point channels
		up0 = np.array( list(itertools.chain.from_iterable(self.point_init_variables) ) )
		len_point_variables = []
		for init, param, init_variables in zip(self.point_init, self.point_params, self.point_init_variables):
			param.update( init(param) )
			len_point_variables.append( len(init_variables) )
		id_start_point_variables = [0]+list(itertools.accumulate(len_point_variables))
		del id_start_point_variables[-1]
		id_end_point_variables = [x+y for (x,y) in zip(id_start_point_variables, len_point_variables)]

		# Init variable
		u0   = np.hstack([ud0, up0])

		# Init event que
		event_que = []
		for i, (param, event, loc) in enumerate(zip(self.point_params, self.point_event, self.point_locations)):
			for etime in param['event_time']:
				if (self.tstart < etime):
					event_que.append( {'time': etime, 'object_id': i, 'func': event, 'loc': loc} )
		event_que = sorted(event_que, key=lambda x: x['time'])
		# print('event_que: ', event_que)

		# Run simulation
		tspan  = [ self.tstart, self.tend ]
		id_que = 0
		final  = 0

		self.t = np.empty(0)                # 0-dim array with size 0
		self.y = np.empty([u0.shape[0], 0]) # 1-dim array with size 0
		while True:
			# Preprocess: Setting of the simulation time up to event-time or to tend.
			if (id_que < len(event_que)) and ( self.tend > event_que[id_que]['time'] ):
				tspan[1] = event_que[id_que]['time']
				final    = 0
			else:
				tspan[1] = self.tend
				final    = 1
			print('tspan: ', tspan)
			print('final: ', final)
			
			# Simulation  (self._f or _ff)
			sol = solve_ivp(self._f, tspan, u0, method='LSODA')
			#sol = self.solve_ivp_numba(tspan, u0, method='LSODA')
			
			# Postprocess: data storage
			self.t = np.hstack([self.t, sol.t])
			self.y = np.hstack([self.y, sol.y])
			if final == 1:
				# End of simulation
				break
			else:
				# Event handling
				t_end = self.t[-1]
				y_end = self.y[:, -1].copy()
				v0, ud0, up0 = np.split(y_end, [self.n_compartments, ud_num])
				loc = event_que[id_que]['loc']
				i   = event_que[id_que]['object_id']
				ids = list(range( id_start_point_variables[i], id_end_point_variables[i] ))
				print('ids: ', ids)
				
				
				# Send the time, memb_pot, and states to event func
				v0[ loc ], up0[ ids ] = event_que[id_que]['func'](t_end, v0[ loc ], up0[ ids ] )
				id_que += 1
				tspan[0] = t_end
				# Set the initial condition.
				u0       = np.hstack([v0, ud0, up0])
		return True


	def _f(self, t, u):

		cm = self.cm
		rm = self.rm
		ri_inv = self.ri_inv
		connector_to_compartments = self.connector_to_compartments
		compartment_to_connectors = self.compartment_to_connectors
		n_compartments = self.n_compartments
		n_connectors   = self.n_connectors
		ud_num = self.ud_num

		# Variables:
		# v : membrane potentials at compartments
		# ud: internal states (distributed mechanisms)
		# up: internal states (point mechanisms)

		v, ud, up = np.split(u, [n_compartments, ud_num])
		ud        = ud.reshape((-1, n_compartments))
		dv        = np.zeros_like(v)
		# Passive channel
		dv += - (v-self.p['El']) / rm # passive property

		# Voltage transfer among compartments
		v_connectors = np.zeros(n_connectors, dtype=float)
		for id_con, id_comps in connector_to_compartments.items() :

			# print('id_comps', id_comps, ' np.sum(ri_inv[id_comps]) ', np.sum(ri_inv[id_comps]) )
			v_connectors[id_con] = np.sum( v[id_comps] * ri_inv[id_comps] ) / np.sum(ri_inv[id_comps])
		for id_comp, id_cons in compartment_to_connectors.items() :
			dv[id_comp]  += \
				self.ri_inv[id_comp] * ( np.sum(v_connectors[id_cons] - v[id_comp]) )

		# Distributed mechanisms
		i_s = 0
		d_ud = np.empty_like(ud)
		for mechanism, num_variables in zip(self.distributed_mechanisms, self.distributed_num_variables):
			
			Im, dd_ud = mechanism(v, dv, ud[i_s:i_s+num_variables, :])
			d_ud[i_s:i_s+num_variables, :] = dd_ud
			dv  += Im
			i_s += num_variables

		# Point mechanisms
		i_s = 0
		I = np.zeros_like(v)
		d_up = np.empty_like(up)
		for mechanism, param, num_variables, loc in \
			zip(self.point_mechanisms, self.point_params, self.point_num_variables, self.point_locations):
			I, d_up[i_s:i_s+num_variables] = mechanism(t, v[loc], up[i_s:i_s+num_variables], param)
			i_s += num_variables
			dv[loc] += I

		# Altogether
		dv   = dv / self.cm
		dudt = np.vstack((dv, d_ud)).reshape(-1)
		dudt = np.hstack([dudt, d_up])

		return dudt



if __name__ == '__main__':
	pass

