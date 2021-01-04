#!/usr/bin/env python3
#
# Strength-duration curves
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import numpy as np
import os
import shared
import sys

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-6'

# Load model
model = shared.model('voigt')
name = model.name()
print(f'Loaded {name}')

# Gather data
spath = 'sd-data'
ks = shared.ko_levels
cs = shared.ko_colors
sd1 = []
sd2 = []
for level in ks:
    level = float(level)
    level_str = str(level).replace('.', '-')

    path = os.path.join(spath, f'sd-{name}-ko-{level_str}-b-1.csv')
    if os.path.isfile(path):
        print(f'Reading {path}')
        sd1.append(myokit.DataLog.load_csv(path))
    else:
        print(f'File not found: {path}')
        sd1.append(None)

    path = os.path.join(spath, f'sd-{name}-ko-{level_str}-b-900.csv')
    if os.path.isfile(path):
        print(f'Reading {path}')
        sd2.append(myokit.DataLog.load_csv(path))
    else:
        print(f'File not found: {path}')
        sd2.append(None)

#
# Create figure
#
fig = plt.figure(figsize=(4.3, 4.5))  # One-column size
fig.subplots_adjust(0.15, 0.100, 0.98, 0.945)
grid = matplotlib.gridspec.GridSpec(2, 1, hspace=0.2)

# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.975, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.004, 0.920, 'A', letter_font)
fig.text(0.004, 0.480, 'B', letter_font)


# Early
ax = fig.add_subplot(grid[0, 0])
ax.set_ylabel('-1x stimulus\namplitude (A/F)')
ax.set_xlim(0, 5)
ax.set_ylim(0, 20)
for k, c, d in zip(ks, cs, sd1):
    if d is not None:
        ax.plot(d['duration'], d['amplitude'], color=c, label=f'{k} mM')
#ax.text(0.96, 0.83, atext, transform=ax.transAxes, ha='right')
ax.legend(loc='upper right')
ax.text(0.10, 1, 'After 1 second')

# Late
ax = fig.add_subplot(grid[1, 0])
ax.set_xlabel('Stimulus duration (ms)')
ax.set_ylabel('-1x stimulus\namplitude (A/F)')
ax.set_xlim(0, 5)
ax.set_ylim(0, 20)
for k, c, d in zip(ks, cs, sd2):
    if d is not None:
        ls = '--' if k == 5.4 else '-'
        ax.plot(d['duration'], d['amplitude'], color=c, ls=ls)
ax.text(0.10, 1, 'After 15 minutes')

#
# Store
#
fname = f'{fname}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')

