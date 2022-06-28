import numpy as np



class PointNMDAR():

	def __init__(self):
		self.p['gNMDAR (default)']        = 2.0
		self.p['ENMDAR (default)']        = 0.0
		self.p['t_act_NMDAR (default)']   = 0.005
		self.p['t_inact_NMDAR (default)'] = 0.05

	def insert_nmdar(self, loc):
		p = {}
		p['gNMDAR']        = self.p['gNMDAR (default)']
		p['ENMDAR']        = self.p['ENMDAR (default)']
		p['t_act_NMDAR']   = self.p['t_act_NMDAR (default)']
		p['t_inact_NMDAR'] = self.p['t_inact_NMDAR (default)']
		p['event_time']    = []
		self.point_params.append(p)
		self.point_locations.append(loc)
		self.point_mechanisms.append(self.nmdar)
		self.point_init.append(self.nmdar_init)
		self.point_event.append(self.nmdar_event)
		self.point_num_variables.append(2)
		self.point_init_variables.append([0, 0])
		return p

	def nmdar_init(self, p):
		t_act_NMDAR   = p['t_act_NMDAR']
		t_inact_NMDAR = p['t_inact_NMDAR']
		t_max         = t_act_NMDAR * t_inact_NMDAR / (t_inact_NMDAR - t_act_NMDAR)
		t_max         *= ( np.log(t_inact_NMDAR) - np.log(t_act_NMDAR))
		c_max         = np.exp(-t_max/t_inact_NMDAR) - np.exp(-t_max/t_act_NMDAR)
		p['t_max']    = t_max
		p['c_max']    = c_max

	def nmdar_event(self, t, v, s):
		s += 1
		return v, s

	def nmdar(self, t, v, s, p):
		gNMDAR        = p['gNMDAR']
		ENMDAR        = p['ENMDAR']
		t_act_NMDAR   = p['t_act_NMDAR']
		t_inact_NMDAR = p['t_inact_NMDAR']
		c_max         = p['c_max']
		ds = np.zeros_like(s)
		ds[0] = -1/t_act_NMDAR   * s[0]
		ds[1] = -1/t_inact_NMDAR * s[1]
		I  = gNMDAR / c_max * ( s[1] - s[0] ) * ( ENMDAR - v )
		return I, ds


class PointCurrent():

	def __init__(self):
		self.p['I (default)']        = 10.0
		self.p['t_start (default)']  = 10.0
		self.p['t_stop (default)']   = 40.0

	def insert_current(self, loc):
		p = {}
		p['I']       = self.p['I (default)']
		p['t_start'] = self.p['t_start (default)']
		p['t_stop']  = self.p['t_stop (default)']
		self.point_params.append(p)
		self.point_locations.append(loc)
		self.point_mechanisms.append(self.current)
		self.point_init.append(self.current_init)
		self.point_event.append(self.current_event)
		self.point_num_variables.append(0)
		self.point_init_variables.append([])
		return p

	def current_init(self, p):
		#return p
		pass

	def current_event(self, t, v, s):
		return v, s

	def current(self, t, v, s, p):

		t_start  = p['t_start']
		t_stop   = p['t_stop']

		if (t < t_stop and t >= t_start):
			I = p['I']
		else:
			I = 0

		ds = np.zeros_like(s)

		return I, ds



class DistributedNa():
	def insert_na(self):
		self.na_init_variables = [0.05, 0.60]
		self.p['gNa']=120.0
		self.p['ENa']= 55.0
		self.distributed_mechanisms.append(self.na)
		self.distributed_num_variables.append( 2)
		self.distributed_init.append(self.na_init)

	def na_init(self, ud0):
		ud0.extend( self.na_init_variables )

	def na(self, v, dv, s):
		gNa= self.p['gNa']
		ENa= self.p['ENa']
		ds = np.zeros_like(s)

		# Na channel
		am = 0.1   * (v+35.0)/(1-np.exp(-(v+35.0)/10.0))
		bm = 4.0   * np.exp(-(v+60.0)/18.0)
		ah = 0.07  * np.exp(-(v+60.0)/20.0)
		bh = 1/(np.exp(-(v+30.0)/10.0)+1.0)

		m = s[0,:]
		h = s[1,:]
		ds[0,:] = am * (1-m) -bm * m
		ds[1,:] = ah * (1-h) -bh * h
		dv += - gNa*m*m*m*h*(v-ENa)

		return ds



class DistributedK():
	def insert_k(self):
		self.k_init_variables = [0.3]
		self.p['gK'] = 36.0
		self.p['EK'] = -72.0
		self.distributed_mechanisms.append(self.k)
		self.distributed_num_variables.append(1)
		self.distributed_init.append(self.k_init)

	def k_init(self, ud0):
		ud0.extend( self.k_init_variables )

	def k(self, v, dv, s):
		gK = self.p['gK']
		EK = self.p['EK']
		ds = np.zeros_like(s)

		# K channel
		n = s[0,:]
		an = 0.01  * (v+50.0)/(1-np.exp(-(v+50.0)/10.0))
		bn = 0.125 * np.exp(-(v+60.0)/80.0)
		ds = an * (1-n) -bn * n
		dv +=-gK*n*n*n*n*(v-EK)

		return ds

