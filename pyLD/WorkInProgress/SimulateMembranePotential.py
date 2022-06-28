# from numbalsoda import lsoda_sig, lsoda
# from numba import njit, cfunc
# import numba as nb
# from numba import njit, jit, float64, int64, f8, i8, b1, void, literal_unroll


import numpy as np
import itertools
from scipy.integrate import solve_ivp
from MembraneMechanisms import DistributedNa, DistributedK, PointCurrent, PointNMDAR



class SimulateMembranePotential(DistributedNa, DistributedK, PointNMDAR, PointCurrent):
	"""Simulate membrane potential of a multicompartment model. This was purely written in Python to simulate a multicompartment model that has a realistic shape. Users can insert distributed and point mechanisms. The distributed mechanisms denotes the distributed channels of Na, K, and Ca. Only one instance for each mechanism can be generated. The point mechanisms denotes the currents via NMDARs and patch pipettes. Multple point mechanisms can be inserted simultaneously.

	Args:
		None

	Returns:
		(pyLD.SimulateMembranePotential): SimulateMembranePotential object that has the follwing instances:
		- init_V (float): Initial membrane potential
		- p (dict): Parameters of pasive properties 'Cm' (XX), 'gL' (XX), and 'EL' (mV).
		- connection (tuple): Connections between compartments (1D array; [id1, id2, id3, ...]).
		- n_compartments (int): Number of compartments.
		- s (numpy[float]): Surface areas of compartments (1D array; [area1, area2, area3, ...]).
		- ra_inv (numpy[float]): Inverse of axial resistances between compartments (1D array; [1?R1, 1/R2, 1/R3, ...]).
		- tstart (float): Start time of simulation.
		- tend (float): End time of simulation.
	"""

	def __init__(self):
		# Passive properties
		self.init_V         = -60.0
		self.p              = {'Cm':1.0, 'gL': 0.3, 'EL': -49.387} # v
		# Compartment
		self.connection = tuple([np.array([i,i+1]) for i in range(9)])
		self.n_compartments = 10
		# Surface area
		self.s    = np.ones(self.n_compartments).astype(np.float64)
		# Axial registance
		self.ra_inv = np.ones(self.n_compartments).astype(np.float64)*0.5
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
		PointNMDAR.__init__(self)
		PointCurrent.__init__(self)

		self.tstart = 0.0
		self.tend   = 50.0

	def _f(self, t, u, s, ra_inv, connection, ud_num ):
		gL = self.p['gL']
		EL = self.p['EL']
		C  = self.p['Cm']
		n_compartments = self.n_compartments

		# Variables
		v, ud, up = np.split(u, [n_compartments, ud_num])
		ud        = ud.reshape((-1, n_compartments))
		dv        = np.zeros_like(v)
		# Passive channel
		dv += - gL*(v-EL)

		# Connection
		# for ii in literal_unroll(connection):
		for ii in connection:
			v_center = np.sum(v[ii]*ra_inv[ii]) / np.sum(ra_inv[ii])
			dv[ii]  += ra_inv[ii]*(v_center - v[ii])

		# Distributed mechanisms
		i_s = 0
		dud = np.empty_like(ud)
		for mechanism, num_variables in zip(self.distributed_mechanisms, self.distributed_num_variables):
			dud[i_s:i_s+num_variables, :] = mechanism(v, dv, ud[i_s:i_s+num_variables, :])
			i_s += num_variables

		# Point mechanisms
		i_s = 0
		I = np.zeros_like(v)
		dup = np.empty_like(up)
		for mechanism, param, num_variables, loc in \
			zip(self.point_mechanisms, self.point_params, self.point_num_variables, self.point_locations):
			I, dup[i_s:i_s+num_variables] = mechanism(t, v[loc], up[i_s:i_s+num_variables], param)
			i_s += num_variables
			dv[loc] += I

		# Altogether
		dudt = np.vstack((dv/C,dud)).reshape(-1)
		dudt = np.hstack([dudt, dup])

		return dudt


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
			init(ud0)
		ud0 = np.tile(np.array(ud0),(self.n_compartments,1)).T.reshape(-1)
		ud_num = ud0.shape[0]

		# Init point channels
		up0 = np.array( list(itertools.chain.from_iterable(self.point_init_variables) ) )
		len_point_variables = []
		for init, param, init_variables in zip(self.point_init, self.point_params, self.point_init_variables):
			init(param)
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

		self.t = np.empty(0)
		self.y = np.empty([u0.shape[0], 0])
		while True:
			# Preprocess
			if (id_que < len(event_que)) and ( self.tend > event_que[id_que]['time'] ):
				tspan[1] = event_que[id_que]['time']
				final    = 0
			else:
				tspan[1] = self.tend
				final    = 1
			# Simulation
			print('tspan: ', tspan)
			print('final: ', final)
			
			sol = solve_ivp(self._f, tspan, u0, method='RK23', args=(self.s, self.ra_inv, self.connection, ud_num))
			self.t = np.hstack([self.t, sol.t])
			self.y = np.hstack([self.y, sol.y])
			# Postprocess
			if final == 1:
				break
			else:
				t_end = self.t[-1]
				y_end = self.y[:, -1].copy()
				v0, ud0, up0 = np.split(y_end, [self.n_compartments, ud_num])
				loc = event_que[id_que]['loc']
				i   = event_que[id_que]['object_id']
				ids = list(range( id_start_point_variables[i], id_end_point_variables[i] ))
				print('ids: ', ids)
				v0[ loc ], up0[ ids ] = event_que[id_que]['func'](t_end, v0[ loc ], up0[ ids ] )
				
				id_que += 1
				tspan[0] = t_end
				u0       = np.hstack([v0, ud0, up0])
		return True


if __name__ == '__main__':
	pass

