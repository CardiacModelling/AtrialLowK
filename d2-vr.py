#!/usr/bin/env python3
#
# Show literature data Vr versus Ko
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import numpy as np
import os
import shared

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-d2-vr'

# Load and adjust literature data
dpath = 'literature-data'
d1 = os.path.join(dpath, 'gelband-1972-fig-6.csv')
d1 = myokit.DataLog.load_csv(d1).npview()
d2 = os.path.join(dpath, 'ten-eick-1979-fig-4.csv')
d2 = myokit.DataLog.load_csv(d2).npview()

#
# Create figure
#
fig = plt.figure(figsize=(9, 5))  # Two-column size
#fig.subplots_adjust(0.15, 0.15, 0.98, 0.98)
fig.subplots_adjust(0.1, 0.1, 0.98, 0.98)
grid = matplotlib.gridspec.GridSpec(1, 1, wspace=0.5)

# Minor ticks
#matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True

ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel('[K$^+$]$_o$ (mM)')
ax.set_ylabel('$V_r$ (mV)')
ax.set_xscale('log')
ax.plot(d1['Ko'], d1['Vr'], 's-', label='Gelband et al. 1972 (n=6)')
ax.plot(d2['Ko'], d2['Vr'], 'x-', label='Ten Eick et al. 1979 (n=2)')
ax.set_xticks(d1['Ko'])
ax.set_xticklabels([int(x) for x in d1['Ko']])
ax.legend(loc='upper left')

#
# Store
#
fname = f'{fname}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')

