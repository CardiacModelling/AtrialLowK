#!/usr/bin/env python3
#
# Supplementary figure: AP and CaT for all models
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-a1-ap-cat'
debug = 'debug' in sys.argv

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare models
models = list(shared.models().values())
for model in models:
    shared.prepare_model(model, protocol, pre_pace=not debug)

# Maximum time to show in plots
tmax = 500

# Create figure
fig = plt.figure(figsize=(9, 9))    # Two-column size
fig.subplots_adjust(0.08, 0.06, 0.98, 0.99, hspace=0.2)
grid = matplotlib.gridspec.GridSpec(3, 1)

ax1 = fig.add_subplot(grid[0:2, 0])
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('V (mV)')

ax2 = fig.add_subplot(grid[2, 0])
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('[Ca]i (uM)')

for model in models:
    label = shared.fancy_name(model)
    vm = model.labelx('membrane_potential')
    cai = model.label('Ca_i')
    if cai is None:
        print(f'No Ca_i marked for {label}.')
    else:
        cai.convert_unit('uM')

    s = myokit.Simulation(model, protocol)
    d = s.run(tmax)
    ax1.plot(d.time(), d[vm], label=label)
    if cai is not None:
        ax2.plot(d.time(), d[cai], label=label)
    else:
        ax2.plot([0], [0], ':')  # Keep color cycle in sync

ax1.legend(loc='upper right')
ax2.legend(loc='upper right')



fig.align_labels()

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
