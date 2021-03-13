#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents
#
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import myokit.lib.plots as mp
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s3-relative-contributions'
debug = 'debug' in sys.argv
shared.splash(fname)

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare models
models = shared.models()
for model in models.values():
    shared.prepare_model(
        model, protocol, fix_cleft_ko=False, pre_pace=not debug)

# Maximum time to show in plots
tmax = 800


def text(ax, x, y, t, c='w'):
    ax.text(x, y, t, color=c, transform=ax.transAxes, fontweight='bold',
            horizontalalignment='right', verticalalignment='center')


# Create figure
fig = plt.figure(figsize=(9, 9))
fig.subplots_adjust(0.075, 0.05, 0.98, 0.97, hspace=0.35, wspace=0.2)
grid = matplotlib.gridspec.GridSpec(3, 3)

#
# Top row: Grandi models
#
# Voigt-Heijman 2013
model = models['voigt']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 0])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)
text(ax, 0.950, 0.865, 'INaK')
text(ax, 0.950, 0.625, 'IK1')
text(ax, 0.950, 0.423, 'ICl,B')
text(ax, 0.950, 0.306, 'INaCa')
text(ax, 0.950, 0.189, 'ICa,B')
text(ax, 0.950, 0.061, 'INa,B')
text(ax, 0.207, 0.542, 'IKur')
text(ax, 0.230, 0.440, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.202, 0.900, -0.098, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.510, 0.640, 'IKr')
ax.arrow(0.410, 0.640, -0.086, -0.030, transform=ax.transAxes, zorder=999)

# Grandi 2011
model = models['grandi']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 1])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)
text(ax, 0.950, 0.865, 'INaK')
text(ax, 0.950, 0.625, 'IK1')
text(ax, 0.950, 0.434, 'ICl,B')
text(ax, 0.950, 0.316, 'INaCa')
text(ax, 0.950, 0.190, 'ICa,B')
text(ax, 0.950, 0.065, 'INa,B')
text(ax, 0.207, 0.542, 'IKur')
text(ax, 0.230, 0.440, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.202, 0.900, -0.098, -0.120, transform=ax.transAxes, zorder=999)


#
# Middle row: Courtemanche models
#
# Courtemanche 1998
model = models['courtemanche']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[1, 0])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)
text(ax, 0.950, 0.950, 'ICa,P')
text(ax, 0.950, 0.838, 'INaK')
text(ax, 0.950, 0.630, 'IK1')
text(ax, 0.950, 0.440, 'INaCa')
text(ax, 0.950, 0.236, 'ICa,B')
text(ax, 0.950, 0.055, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.210, 0.450, 'ICaL')
text(ax, 0.270, 0.840, 'Ito')
ax.arrow(0.172, 0.840, -0.088, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.440, 0.660, 'IKs')
ax.arrow(0.332, 0.660, -0.093, -0.060, transform=ax.transAxes, zorder=999)
text(ax, 0.420, 0.580, 'IKr')
ax.arrow(0.320, 0.580, -0.086, -0.030, transform=ax.transAxes, zorder=999)

# Ni 2017
model = models['ni']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[1, 1])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)
text(ax, 0.950, 0.950, 'ICa,P')
text(ax, 0.950, 0.838, 'INaK')
text(ax, 0.950, 0.635, 'IK1')
text(ax, 0.950, 0.420, 'INaCa')
text(ax, 0.950, 0.230, 'ICa,B')
text(ax, 0.950, 0.055, 'INa,B')
text(ax, 0.215, 0.550, 'IKur')
text(ax, 0.225, 0.440, 'ICaL')
text(ax, 0.280, 0.850, 'Ito')
ax.arrow(0.182, 0.850, -0.088, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.500, 0.720, 'IKs')
ax.arrow(0.392, 0.720, -0.093, -0.060, transform=ax.transAxes, zorder=999)
text(ax, 0.480, 0.640, 'IKr')
ax.arrow(0.380, 0.640, -0.086, -0.030, transform=ax.transAxes, zorder=999)

#
# Bottom row: Nygren models
#
# Nygren 1998
model = models['nygren']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 0])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

text(ax, 0.950, 0.870, 'INaK')
text(ax, 0.950, 0.640, 'IK1')
text(ax, 0.950, 0.420, 'INaCa')
text(ax, 0.950, 0.230, 'ICa,B')
text(ax, 0.950, 0.060, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.220, 0.450, 'ICaL')
text(ax, 0.300, 0.880, 'Ito')
ax.arrow(0.20, 0.88, -0.11, -0.08, transform=ax.transAxes, zorder=999)
text(ax, 0.400, 0.06, 'INa')
ax.arrow(0.28, 0.06, -0.09, -0.03, transform=ax.transAxes, zorder=999)

# Maleckar 2009
model = models['maleckar']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 1])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)
text(ax, 0.950, 0.870, 'INaK')
text(ax, 0.950, 0.620, 'IK1')
text(ax, 0.950, 0.445, 'INaCa')
text(ax, 0.950, 0.320, 'ICa,B')
text(ax, 0.950, 0.125, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.220, 0.450, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.20, 0.90, -0.11, -0.08, transform=ax.transAxes, zorder=999)
text(ax, 0.380, 0.06, 'INa')
ax.arrow(0.26, 0.06, -0.09, -0.03, transform=ax.transAxes, zorder=999)

# Koivumaki 2011
model = models['koivumaki']
currents, colours = shared.currents(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 2])
ax.set_title(shared.fancy_name(model))
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)
text(ax, 0.950, 0.880, 'INaK')
text(ax, 0.950, 0.640, 'IK1')
text(ax, 0.950, 0.440, 'INaCa')
text(ax, 0.950, 0.260, 'ICa,B')
text(ax, 0.950, 0.080, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.215, 0.450, 'ICaL')
text(ax, 0.300, 0.880, 'Ito')
ax.arrow(0.20, 0.88, -0.11, -0.08, transform=ax.transAxes, zorder=999)
text(ax, 0.360, 0.053, 'INa')
ax.arrow(0.24, 0.053, -0.09, -0.03, transform=ax.transAxes, zorder=999)

#
# Legend
#
ax = fig.add_subplot(grid[0, 2])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
lines = []
for i, current in enumerate(shared.current_names):
    lines.append(matplotlib.lines.Line2D([0], [0], color=shared.cmap(i), lw=5))
labels = [x.replace('_', '') for x in shared.current_names]
#ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)


# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
