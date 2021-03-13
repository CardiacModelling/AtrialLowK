#!/usr/bin/env python3
#
# Rectification as a function of voltage and Ko
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
fname = 'figure-s6-ik1-x'
shared.splash(fname)

many = True


def x_voigt(V, K=5.4):
    """ Calculate ``x_infty`` for the Voigt-Heijman model. """
    Ki = 120    # mM
    RTF = 8314 * 310 / 96485    # mV
    E = RTF * np.log(K / Ki)
    Na = 9.4
    a = 0.1 + 0.9 / (1 + (Na / 7)**2)
    a /= 1 + np.exp(0.2385 * (V - E - 59.215))
    b = 0.49124 * np.exp(0.08032 * (V - E + 5.476))
    b += np.exp(0.06175 * (V - E - 594.31))
    b /= 1 + np.exp(-0.5143 * (V - E + 4.753))
    return a / (a + b)


def x_grandi(V, K=5.4):
    """ Calculate ``x_infty`` for the Grandi / LR1 model. """
    Ki = 120    # mM
    RTF = 8314 * 310 / 96485    # mV
    E = RTF * np.log(K / Ki)
    a = 1.02
    a /= 1 + np.exp(0.2385 * (V - E - 59.215))
    b = 0.49124 * np.exp(0.08032 * (V - E + 5.476))
    b += np.exp(0.06175 * (V - E - 594.31))
    b /= 1 + np.exp(-0.5143 * (V - E + 4.753))
    return a / (a + b)


def x_courtemanche(V, K=5.4):
    """ Calculate ``x_infty`` for the Courtemanche model. """
    return 1 / (1 + np.exp(0.07 * (V + 80)))


def x_nygren(V, K=5.4):
    """ Calculate ``x_infty`` for the Nygren model. """
    Ki = 129.435
    RTF = 8314 * 306.15 / 96487     # mV
    E = RTF * np.log(K / Ki)
    return 1 / (1 + np.exp(1.5 * (V - E + 3.6) / RTF))


def plot(ax, f, title, default_only=False):
    ax.axhline(0, ls='--', color='#cccccc')
    ax.axhline(1, ls='--', color='#cccccc')
    ax.set_xlim(-130, 0)
    ax.set_xlabel('V (mV)')
    ax.set_ylabel('x (-)')
    if default_only:
        k = shared.ko_levels[3]
        c = shared.ko_colors[3]
        ax.plot(V, f(V, k), c, label=f'{k} mM')
    else:
        if many:
            for k in np.linspace(0, 10, 40):
                ax.plot(V, f(V, k), color='k', alpha=0.1)
        for k, c in zip(shared.ko_levels, shared.ko_colors):
            ax.plot(V, f(V, k), c, label=f'{k} mM')
    ax.legend(loc='center right')
    ax.text(1, 0.85, title, horizontalalignment='right',
            transform=ax.transAxes)


#
# Create figure
#
fig = plt.figure(figsize=(9, 5))  # Two-column size
fig.subplots_adjust(0.06, 0.09, 0.99, 0.99)
grid = matplotlib.gridspec.GridSpec(2, 2, wspace=0.32, hspace=0.30)
V = np.linspace(-140, 20, 200)
ax = fig.add_subplot(grid[0, 0])
plot(ax, x_voigt, 'Voigt-Heijman et al.')
ax = fig.add_subplot(grid[0, 1])
plot(ax, x_nygren, 'Nygren et al.')
ax = fig.add_subplot(grid[1, 0])
plot(ax, x_grandi, 'Grandi-Pandit-Voigt et al. (LR1)')
ax = fig.add_subplot(grid[1, 1])
plot(ax, x_courtemanche, 'Courtemanche et al.', True)


# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
