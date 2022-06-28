
from SimulateMembranePotential import SimulateMembranePotential

m = SimulateMembranePotential()
m.insert_na()
m.insert_k()
p1 = m.insert_nmdar(0)
p1['event_time'] = [10,20,30]
p2 = m.insert_nmdar(1)
p2['event_time'] = [5,20,45]
# p = m.insert_current(0)
m.run()

import matplotlib.pyplot as plt

for i in range(m.n_compartments):
	plt.plot(m.t,m.y[i,:],'-', label=str(i))

plt.legend()
plt.show()


