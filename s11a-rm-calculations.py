#!/usr/bin/env python3
#
# Rm calculcations with different steps
#
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s11a-rm'
debug = 'debug' in sys.argv
shared.splash(fname)

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare models
models = shared.models()
for name, model in models.items():
    shared.prepare_model(
        model, protocol, fix_cleft_ko=False, pre_pace=not debug)

# Maximum time to show in plots
tmax = 800

# Create figure
fig = plt.figure(figsize=(9, 10))
fig.subplots_adjust(0.07, 0.045, 0.92, 0.99, hspace=0.095, wspace=0.09)
grid = matplotlib.gridspec.GridSpec(len(models), 4)


def rm_plot(fig, i, j, model, protocol, simulation, t, v, dt, dv):
    """Plot AP and Rm"""
    simulation.reset()
    d = simulation.run(tmax, log=[t, v])

    ax = fig.add_subplot(grid[i, j])
    if i == len(models) - 1:
        ax.set_xlabel('Time (ms)')
    else:
        ax.set_xticklabels('')
    if j == 0:
        ax.set_ylabel('V (mV)')
    else:
        ax.set_yticklabels('')
    ax.plot(d.time(), d[v], color='k')

    ax = ax.twinx()
    if j == 3:
        ax.set_ylabel('Rm (MOhm)')
        ax.yaxis.set_label_position('right')
    else:
        ax.set_yticklabels('')
    ax.set_ylim(0, 2000)
    ax.text(0.95, 0.85, f'dt={dt} dv={dv}',
            horizontalalignment='right',
            fontdict={'size': 8},
            transform=ax.transAxes)
    if j == 0:
        ax.text(0.95, 0.70, model.name(),
                horizontalalignment='right',
                fontdict={'size': 7},
                transform=ax.transAxes)

    if not debug:
        ts, rs = shared.rm(model, protocol, dt, dv, tmax, simulation)
        ax.plot(ts, rs, color='tab:red', label=f'{dt}ms {dv}mV')


# Plot for all
for i, model in enumerate(models.values()):
    print('Adding plots for ' + model.name())

    t = model.time()
    v = model.label('membrane_potential')

    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    rm_plot(fig, i, 0, model, protocol, s, t, v, 5, 5)
    rm_plot(fig, i, 1, model, protocol, s, t, v, 10, 5)
    rm_plot(fig, i, 2, model, protocol, s, t, v, 20, 2)
    rm_plot(fig, i, 3, model, protocol, s, t, v, 0.1, 10)

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
