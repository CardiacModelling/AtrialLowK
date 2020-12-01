#!/usr/bin/env python3
#
# Model figure 1: Short-term changes due [K+]o changes
#
# Panel A:
#   AP versus time. 1Hz pacing, only 800ms shown.
#   Pre-paced with baseline Ko=5.4, then switched to test level, then paced for
#   1 beat before plotting
# Panel B:
#   IV relationship for IK1 at each level
# Panel C:
#   Rm at each level
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from  matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-m1-ko'
debug = 'debug' in sys.argv

# Load model
model = shared.model()
name = model.name()

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 800

# Time to pre-pace after changing [K]o level
pre_time = 1000

# Prepare model
shared.prepare_model(model, protocol, fix_cleft_ko=True, pre_pace=not debug)

# Get some variables
v = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Perform simulations
ks = shared.ko_levels if not debug else [5.4]
cs = shared.ko_colors if not debug else [None]
ds = [] # Logs with an AP in
rm = [] # Vectors time, Rm
vr = [] # Vrest, after 1 second

for k in ks:
    print(f'Simulating {k}mM')

    # Create new simulation (otherwise pre() can't be used!)
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    s.set_constant(ko, k)
    if pre_time:
        s.pre(pre_time)
    vr.append(s.state()[v.indice()])
    d = s.run(tmax).npview()
    ds.append(d)

    if not debug:
        dt = 0.1
        dv = 10
        rm.append(shared.rm(model, protocol, dt, dv, tmax, simulation=s))
    del(s)

# Calculate IV curves
fixed = model.clone()
fki = fixed.labelx('K_i')
fko = fixed.labelx('K_o')
fik1 = fixed.labelx('I_K1')
fixed_states = list(fixed.states())
for state in fixed_states:
    if state.label() != 'membrane_potential':
        shared.demote(state)
iv = []
vs = np.linspace(-120, 20, 141)
ivr = []    # I at Vrest
for i, k in enumerate(ks):
    fko.set_rhs(k)
    f = fik1.pyfunc()
    iv.append((vs, f(vs)))
    ivr.append(f(vr[i]))

#
# Create figure
#
fig = plt.figure(figsize=(9, 2.7))  # Two-column size
fig.subplots_adjust(0.08, 0.17, 0.95, 0.96)
grid = GridSpec(1, 3, wspace=0.4)

# Add model name
fig.patches.append(Rectangle((0.965, 0), 0.035, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.975, 0.5, name, name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.003, 0.91, 'A', letter_font)
fig.text(0.320, 0.91, 'B', letter_font)
fig.text(0.633, 0.91, 'C', letter_font)


# A: Change in AP waveform
ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel('Time (ms)')
ax.set_ylabel('V (mV)')
for k, d, c in zip(ks, ds, cs):
    t = pre_time + d.time()
    ax.plot(t, d[v], label=str(k) + ' mM', color=c)
ax.legend(loc='upper right')


# B: IK1 dependence on [K]o
ax = fig.add_subplot(grid[0, 1])
ax.set_xlabel('V (mV)')
ax.set_ylabel('IK1 (A/F)')
ax.axhline(0, color='#bbbbbb')
for xy, c in zip(iv, cs):
    ax.plot(xy[0], xy[1], color=c)
ax.set_xlim(-112, -12)
if 'grandi-20' in name:
    ax.set_ylim(-0.3, 0.3)
else:
    ax.set_ylim(-1.5, 1.0)

for x, y, c in zip(vr, ivr, cs):
    ax.plot(x, y, 'o', markersize=10, fillstyle='none', color=c)

# C: Change in input resistance
ax = fig.add_subplot(grid[0, 2])
ax.set_xlabel('Time (ms)')
ax.set_ylabel('R (MOhm)')
ax.set_ylim(-400, 16000)
if name == 'grandi-2011':
    axi = ax.inset_axes([0.37, 0.37, 0.6, 0.6])
    axi.set_ylim(0, 950)
if not debug:
    for tr, c in zip(rm, cs):
        ts, rs = tr
        t = (pre_time + ts)
        ax.plot(t, rs, color=c)
        if name == 'grandi-2011':
            axi.plot(t, rs, color=c)


# Show / store
path = os.path.join(fpath, fname + '-' + name)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
