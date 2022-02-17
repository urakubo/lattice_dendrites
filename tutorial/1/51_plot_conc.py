import os, sys, glob
import numpy as np
from pyLD import *
import matplotlib.pyplot as plt

m  = 'YFP'
plot_filename  = 'imgs/profile_photobleach.png'
simulation_dir = 'results_photobleach'
t_offset       = -4

'''
m  = 'Ca'
plot_filename  = 'imgs/profile_Ca.png'
simulation_dir = 'results_Ca_dynamics'
t_offset       = -2
'''

label_file     = 'models/labels_ball_and_stick.h5'
lm_files       = sorted( glob.glob(os.path.join(simulation_dir, '*.lm')) )
conc_files     = [f[:-3]+'.h5' for f in lm_files]


print('\nConnect concs at labeled regions.')
s = GetLabeledConcs()
s.load_label_volume(label_file)
for (lm_file, conc_file) in zip(lm_files, conc_files):
    s.lm_file = lm_file
    s.exec()
    s.save(conc_file)

c = ConnectLabeledConcs(conc_files)

print('\nConnect total concs.')
domain_id = 1 # cytosol
t = ConnectTotalConcs(lm_files, 'cytosol')

print('\nPlot figure')
fig = plt.figure(figsize=(6,4))
ax  = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.plot(c.timepoints+t_offset,\
        c.get_concs( species=m, label_ids=1 ),\
        label = m+' in spine' )

ax.plot(t.timepoints+t_offset,\
        t.get_concs( species=m ),\
        label = m+' in total' )

# ax.set_ylim([0, 1.5])
ax.set_xlim([t_offset,np.max(t.timepoints+t_offset)])
plt.xlabel('Time (s)')
plt.ylabel('Conc (uM)')
hans, labs = ax.get_legend_handles_labels()
ax.legend(handles=hans, labels=labs, frameon=False)
plt.savefig(plot_filename, dpi=150)
plt.show()

