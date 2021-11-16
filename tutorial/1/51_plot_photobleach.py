import os, sys
import numpy as np
from pyLD import *
import matplotlib.pyplot as plt

m  = 'YFP'
input_label_file = 'models/labels_ball_and_stick.h5'
output_plot_file = 'imgs/timecourse_photobleach.png'
lm_files         = ['results_photobleach/0000.lm'  ,'results_photobleach/0001.lm']
conc_files       = ['results_photobleach/c0000.h5' ,'results_photobleach/c0001.h5']

# Save and connect concs of labeled spines
s = GetLabeledConcs()
s.load_label_volume(input_label_file)
for (lm_file, conc_file) in zip(lm_files, conc_files):
    s.lm_file = lm_file
    s.exec()
    s.save(conc_file)

c = ConnectLabeledConcs(conc_files)


# Connect total concs
domain_ids = 1
t = ConnectTotalConcs(lm_files, domain_ids)


# Plot a figure
fig = plt.figure(figsize=(6,4))
ax  = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

t_offset = 4

ax.plot(c.timepoints-t_offset,\
        c.get_concs( species=m, label_ids=1 ),\
        label = m + ' in spine' )

ax.plot(t.timepoints-t_offset,\
        t.get_concs( species=m ),\
        label = m+' in total' )

# ax.set_ylim([0, 1.5])
ax.set_xlim([-t_offset,np.max(t.timepoints-t_offset)])
plt.xlabel('Time (s)')
plt.ylabel('Conc (uM)')
hans, labs = ax.get_legend_handles_labels()
ax.legend(handles=hans, labels=labs, frameon=False)
plt.savefig(output_plot_file, dpi=150)
plt.show()