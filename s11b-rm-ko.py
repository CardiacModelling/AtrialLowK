#!/usr/bin/env python3
#
# Rm calculcations with different step sizes and durations
#
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s11b-rm-ko'
debug = 'debug' in sys.argv
shared.splash(fname)

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare model
model = shared.model('voigt')
shared.prepare_model(model, protocol, pre_pace=not debug)

# Maximum time to show in plots
tmax = 800

# Time to pre-pace
pre_time = 1000

# Create figure
fig = plt.figure(figsize=(9, 11.6))
fig.subplots_adjust(0.07, 0.045, 0.92, 0.99, hspace=0.095, wspace=0.09)
grid = matplotlib.gridspec.GridSpec(7, len(shared.ko_levels))


def rm_plot(fig, i, j, model, protocol, simulation, color, t, v, dt, dv):
    """Plot AP and Rm"""
    simulation.reset()
    d = simulation.run(tmax, log=[t, v]).npview()

    ax = fig.add_subplot(grid[i, j])
    if i == 6:
        ax.set_xlabel('Time (ms)')
    else:
        ax.set_xticklabels('')
    if j == 0:
        ax.set_ylabel('V (mV)')
    else:
        ax.set_yticklabels('')
    ax.plot(pre_time + d.time(), d[v], color=color)

    ax = ax.twinx()
    if j == 3:
        ax.set_ylabel('Rm (MOhm)')
        ax.yaxis.set_label_position('right')
    else:
        ax.set_yticklabels('')
    ax.set_ylim(0, 3000)
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
        ax.plot(pre_time + ts, rs, color='tab:green', ls=':')
        rs[rs < 0] = float('nan')
        ax.plot(pre_time + ts, rs, color='tab:green', label=f'{dt}ms {dv}mV')


# Plot for all
for i, (k, c) in enumerate(zip(shared.ko_levels, shared.ko_colors)):
    print(f'Adding plots for {k}mM')

    t = model.time()
    v = model.label('membrane_potential')
    model.labelx('K_o').set_rhs(k)

    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    if pre_time:
        s.pre(pre_time)
    rm_plot(fig, 0, i, model, protocol, s, c, t, v, 5, 5)
    rm_plot(fig, 1, i, model, protocol, s, c, t, v, 10, 5)
    rm_plot(fig, 2, i, model, protocol, s, c, t, v, 20, 2)
    rm_plot(fig, 3, i, model, protocol, s, c, t, v, 10, 10)
    rm_plot(fig, 4, i, model, protocol, s, c, t, v, 0.1, 2)
    rm_plot(fig, 5, i, model, protocol, s, c, t, v, 0.1, 5)
    rm_plot(fig, 6, i, model, protocol, s, c, t, v, 0.1, 10)

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
