#!/usr/bin/env python3
#
# dV/dt max versus Ko in single cells
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.gridspec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s9-dvdtmax'
shared.splash(fname)

# Potassium levels
ks = shared.ko_levels + [7, 8, 10, 12]

# Tick labels
tl = [str(int(x)) if int(x) == x else str(x) for x in ks]
tl[1] = None

# Create protocol
cl = 1000
protocol = shared.default_protocol()

# Gather data
dvdt1 = {}
dvdt2 = {}
for model in shared.models().values():

    # Get some variables
    shared.prepare_model(model, protocol)
    v = model.labelx('membrane_potential')
    t = model.time()
    dv = 'dot(' + v.qname() + ')'
    ko = model.labelx('K_o')

    # Perform simulations
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-10, 1e-10)
    dvs1 = []
    dvs2 = []
    for k in ks:
        s.set_constant(ko, k)

        dvdt = np.nan
        s.reset()
        s.set_state(shared.state_fast(model, protocol, k))
        d = s.run(cl, log=[t, v, dv]).npview()
        vmin, vmax = np.min(d[v]), np.max(d[v])
        if vmax > -20 and vmin < -50:
            v90 = 0.9 * vmin + 0.1 * vmax
            apd = d.apd(v, v90)['duration']
            if len(apd) > 0 and apd[0] > 100:
                dvdt = np.max(d[dv])
        dvs1.append(dvdt)

        dvdt = np.nan
        s.reset()
        s.set_state(shared.state_slow(model, protocol, k))
        d = s.run(cl, log=[t, v, dv]).npview()
        vmin, vmax = np.min(d[v]), np.max(d[v])
        if vmax > -20 and vmin < -50:
            v90 = 0.9 * vmin + 0.1 * vmax
            apd = d.apd(v, v90)['duration']
            if len(apd) > 0 and apd[0] > 50:
                dvdt = np.max(d[dv])
        dvs2.append(dvdt)

    dvdt1[model.name()] = dvs1
    dvdt2[model.name()] = dvs2


#
# Create figure
#
fig = plt.figure(figsize=(9, 8))  # Two-column size
fig.subplots_adjust(0.07, 0.065, 0.99, 0.950)
grid = matplotlib.gridspec.GridSpec(3, 2, wspace=0.32, hspace=0.30)

yl = [0, 451]

cmap = matplotlib.cm.get_cmap('tab10')
cs = [cmap(i) for i in range(len(dvdt1))]


def add(row, models, colors):
    ax = fig.add_subplot(grid[row, 0])
    ax.set_xlabel('$[K^+]_o$ (mM)')
    ax.set_ylabel('dV/dt max (V/s)')
    ax.set_xlim(ks[0] - 0.5, ks[-1] + 0.5)
    ax.set_ylim(*yl)
    ax.set_xticks(ks)
    ax.set_xticklabels(tl)
    for m, c in zip(models, colors):
        ax.plot(ks, dvdt1[m.name()], 's-', color=c, label=shared.fancy_name(m))
    #ax.legend()
    if row == 0:
        ax.text(7, 500, 'After 1 second', horizontalalignment='center')

    ax = fig.add_subplot(grid[row, 1])
    ax.set_xlabel('$[K^+]_o$ (mM)')
    ax.set_ylabel('dV/dt max (V/s)')
    ax.set_xlim(ks[0] - 0.5, ks[-1] + 0.5)
    ax.set_ylim(*yl)
    ax.set_xticks(ks)
    ax.set_xticklabels(tl)
    for m, c in zip(models, colors):
        ax.plot(ks, dvdt2[m.name()], 's-', color=c, label=shared.fancy_name(m))
    ax.legend(loc='upper right')
    if row == 0:
        ax.text(7, 5050, 'After 15 minutes', horizontalalignment='center')

models = list(shared.models().values())
add(0, models[:2], cs[:2])
add(1, models[2:4], cs[2:4])
add(2, models[4:], cs[4:])

fig.align_labels()


# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')

