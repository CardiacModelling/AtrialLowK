#!/usr/bin/env python3
#
# Three Ko effects on IK1
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
fname = 'figure-s5-ik1-ko'
shared.splash(fname)

model = None
if len(sys.argv) == 2:
    model = sys.argv[1].strip().lower()
if model not in ['voigt', 'nygren']:
    print(f'Syntax: {sys.argv[0]} model')
    print(f'  where model is voigt or nygren')
    sys.exit(1)


def ik1_voigt(V, K1=5.4, K2=5.4, K3=5.4):
    """
    Plot IK1 for a certain set of Ko levels, using the equations from the
    Voigt et al. model.

    K1: Ko in driving term EK
    K2: Ko in square root dependence
    K3: Ko in kinetics
    """
    Ki = 120    # mM
    RTF = 8314 * 310 / 96485    # mV
    E1 = RTF * np.log(K1 / Ki)
    g = 2.1 * 0.0525 * np.sqrt(K2 / 5.4)
    E3 = RTF * np.log(K3 / Ki)
    Na = 9.4
    a = 0.1 + 0.9 / (1 + (Na / 7)**2)
    a /= 1 + np.exp(0.2385 * (V - E3 - 59.215))
    b = 0.49124 * np.exp(0.08032 * (V - E3 + 5.476))
    b += np.exp(0.06175 * (V - E3 - 594.31))
    b /= 1 + np.exp(-0.5143 * (V - E3 + 4.753))
    return g * a / (a + b) * (V - E1)      # A/F


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


if model == 'voigt':
    ik1 = ik1_voigt
else:
    ik1 = ik1_nygren

#
# Create figure
#
fig = plt.figure(figsize=(9, 9))  # Two-column size
fig.subplots_adjust(0.09, 0.07, 0.98, 0.95)
grid = matplotlib.gridspec.GridSpec(2, 2, wspace=0.32, hspace=0.30)
ws = 0.50
V = np.linspace(-140, 0, 200)

fig.patches.append(matplotlib.patches.Rectangle(
    (0.510, 0.010), 0.480, 0.980, facecolor='#ffefa1',
    fill=True, alpha=0.5,
    transform=fig.transFigure, zorder=-1))
fig.patches.append(matplotlib.patches.Rectangle(
    (0.510, 0.010), 0.480, 0.980, edgecolor='#000000', fill=False,
    transform=fig.transFigure, zorder=-1))
fig.patches.append(matplotlib.patches.Rectangle(
    (0.005, 0.005), 0.990, 0.495, facecolor='#d7ebff', edgecolor='#000000',
    fill=True, alpha=1,
    transform=fig.transFigure, zorder=-2))

fig.text(
    0.25, 0.965, 'Kinetics effects only',
    horizontalalignment='center', verticalalignment='center')
fig.text(
    0.75, 0.965, 'Conductance term and kinetics',
    horizontalalignment='center', verticalalignment='center')
fig.text(
    0.25, 0.47, 'Driving term and kinetics',
    horizontalalignment='center', verticalalignment='center')
fig.text(
    0.75, 0.47, 'All three effects',
    horizontalalignment='center', verticalalignment='center')

def plot(ax, K1, K2, K3):
    ax.axhline(0, color='#bbbbbb')
    ax.set_xlim(-110, -20)
    ax.set_ylim(-0.5, 0.7)
    ax.set_xlabel('V (mV)')
    ax.set_ylabel('I (A/F)')

    d = 5.4
    for k, c in zip(shared.ko_levels, shared.ko_colors):
        I = ik1(V, k if K1 else d, k if K2 else d, k if K3 else d)
        ax.plot(V, I, color=c, label=f'{k} mM')
    txt = 'E' if K1 else r'$\varnothing$'
    txt += ', '
    txt += '$\sqrt{\,}$' if K2 else r'$\varnothing$'
    txt += ', '
    txt += '$x$' if K3 else r'$\varnothing$'
    ax.text(0.4, 0.1, txt, transform=ax.transAxes)


# Top-left: No driving term or conductance
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=grid[0, 0], wspace=ws)
ax = fig.add_subplot(sg[0, 0])
plot(ax, False, False, False)
ax = fig.add_subplot(sg[0, 1])
plot(ax, False, False, True)
ax.legend(loc=(0.3, 0.68))

# Top-right: No driving term, but with conductance
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=grid[0, 1], wspace=ws)
ax = fig.add_subplot(sg[0, 0])
plot(ax, False, True, False)
ax = fig.add_subplot(sg[0, 1])
plot(ax, False, True, True)

# Bottom-left: Driving term, no conductance
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=grid[1, 0], wspace=ws)
ax = fig.add_subplot(sg[0, 0])
plot(ax, True, False, False)
ax = fig.add_subplot(sg[0, 1])
plot(ax, True, False, True)

# Bottom-right: Driving term and conductance
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=grid[1, 1], wspace=ws)
ax = fig.add_subplot(sg[0, 0])
plot(ax, True, True, False)
ax = fig.add_subplot(sg[0, 1])
plot(ax, True, True, True)


# Show / store
path = os.path.join(fpath, fname + '-' + model)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
