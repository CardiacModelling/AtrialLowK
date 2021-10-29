#!/usr/bin/env python3
#
# IK1 rectification
#
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import matplotlib.patches
import shared
import sys

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s5c-ik1'
shared.splash(fname)


def ik1_nygren(V, K1=5.4, K2=5.4, K3=5.4):
    """
    Plot IK1 for a certain set of Ko levels, using the equations from the
    Voigt et al. model.

    K1: Ko in driving term EK
    K2: Ko in square root dependence
    K3: Ko in kinetics
    """
    Ki = 129.435
    RTF = 8314 * 306.15 / 96487     # mV
    E1 = RTF * np.log(K1 / Ki)
    g = 0.03 * K2**0.4457
    E3 = RTF * np.log(K3 / Ki)
    x = 1 / (1 + np.exp(1.5 * (V - E3 + 3.6) / RTF))
    return  g * x * (V - E1)


ik1 = ik1_nygren

#
# Create figure
#
fig = plt.figure(figsize=(9, 4))  # Two-column size
fig.subplots_adjust(0.08, 0.12, 0.98, 0.98)
grid = matplotlib.gridspec.GridSpec(1, 3, wspace=0.5, hspace=0.30)
ws = 0.50
vlim = [-140, 20]
V = np.linspace(vlim[0], vlim[1], 200)


# One: only driving term
ax = fig.add_subplot(grid[0, 0])
ax.axhline(0, color='#999999', lw=1)
ax.set_xlim(*vlim)
ax.set_ylim(-0.3, 0.62)
ax.set_xlabel('$V$ (mV)')
ax.set_ylabel('$I_{K1}$ (A/F)')

d = 5.4
I = ik1(V)
ax.plot(V, I, color='tab:blue')

# Two: add in kinetics
ax = fig.add_subplot(grid[0, 1])

# Three: add in conductance
ax = fig.add_subplot(grid[0, 2])


# Show / store
path = os.path.join(fpath, fname + '-nygren')
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
