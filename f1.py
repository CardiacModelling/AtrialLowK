#!/usr/bin/env python3
#
# Figure 1: Comparison of model with and without Ko-dependence
#
import matplotlib
import matplotlib.gridspec
import matplotlib.patches
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import shared
import sys

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-1'
debug = 'debug' in sys.argv

# Load models
model1 = shared.model('voigt')
model2 = shared.model('courtemanche')

# Times for instantaneous IV
times1 = [20, 100, 200, 320]
times2 = [20, 90, 175, 275]
boxes1 = [(-80, 24), (-30, 22), (-48, 16), (-71, 25)]
boxes2 = [(-105, 38), (-29, 21), (-62, 46), (-100, 53)]

colors = ['#555555', 'tab:blue', 'tab:orange', 'tab:green']
#colors = ['#555555', 'tab:red', 'tab:purple', 'tab:blue']
alpha = 0.3
black = '#444444'

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 800

# Time to pre-pace after changing [K]o level
pre_time = 1000

# Prepare models
shared.prepare_model(model1, protocol, fix_cleft_ko=True, pre_pace=not debug)
shared.prepare_model(model2, protocol, fix_cleft_ko=True, pre_pace=not debug)

# External potassium levels and colours
ks = shared.ko_levels + [8]
cs = shared.ko_colors + ['#777777']
ls = ['-', '-', '-', '-', '--']
if debug:
    ks, cs, ls = [ks[3]], [cs[3]], [ls[3]]

#
# Create figure
#
#fig = plt.figure(figsize=(9, 9))  # Two-column size
#fig.subplots_adjust(0.09, 0.05, 0.98, 0.96)
#grid = matplotlib.gridspec.GridSpec(4, 2, wspace=0.3, hspace=0.4)
fig = plt.figure(figsize=(9, 7.5))  # Two-column size
fig.subplots_adjust(0.08, 0.065, 0.99, 0.95)
grid = matplotlib.gridspec.GridSpec(3, 2, wspace=0.2, hspace=0.3)


# Add model names
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.295, 0.98, shared.fancy_name(model2), name_font)
fig.text(0.795, 0.98, shared.fancy_name(model1), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.003, 0.940, 'A', letter_font)
fig.text(0.003, 0.610, 'B', letter_font)
fig.text(0.003, 0.290, 'C', letter_font)

#
# A. Plot APs
#
def ap(ax, model, protocol, levels):
    # Get some variables
    vm = model.labelx('membrane_potential')
    ko = model.labelx('K_o')

    # Calculate APs, fetch resting potentials
    ds = []
    vr = [] # Vrest, after 1 second
    for k in levels:
        print(f'Simulating {k}mM with {model.name()}')

        # Create new simulation (otherwise pre() can't be used!)
        s = myokit.Simulation(model, protocol)
        s.set_tolerance(1e-8, 1e-8)
        s.set_constant(ko, k)
        if pre_time:
            s.run(pre_time, log=[])
        vr.append(s.state()[vm.indice()])
        d = s.run(tmax).npview()
        ds.append(d)
        del(s)

    # Plot APs
    for k, c, l, d, in zip(ks, cs, ls, ds):
        ax.plot(d.time(), d[vm], label=str(k) + ' mM', color=c, ls=l)

    # Return resting potentials
    return vr

# Plot APs, store Vrest
ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel('Time (ms)')
ax.set_ylabel('V (mV)')
ax.set_ylim(-105, 34)
ax.set_yticks([-100, -80, -60, -40, -20, 0, 20, 40])
points = {5.4: (times2, colors, 's'), 3.2: (times2, colors, 'x')}
vr2 = ap(ax, model2, protocol, ks)
ax.legend(loc='upper right')

ax = fig.add_subplot(grid[0, 1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('V (mV)')
ax.set_ylim(-81, 15)
vr1 = ap(ax, model1, protocol, ks)

#
# B. Calculate and plot IV curves for IK1
#
def ik1_iv(ax, model, protocol, levels, vrests):
    print(f'Calculating IK1 IV for {model.name()}')

    # Calculate from model
    fixed = model.clone()
    fki = fixed.labelx('K_i')
    fko = fixed.labelx('K_o')
    fik = fixed.labelx('I_K1')
    fixed_states = list(fixed.states())
    for state in fixed_states:
        if state.label() != 'membrane_potential':
            shared.demote(state)
    iv = []
    irs = []    # I at Vrest
    vs = np.linspace(-120, 20, 141)
    for k, vrest in zip(levels, vrests):
        fko.set_rhs(k)
        f = fik.pyfunc()
        iv.append((vs, f(vs)))
        irs.append(f(vrest))

    # Plot
    ax.axhline(0, color='#bbbbbb')
    for xy, c, l in zip(iv, cs, ls):
        ax.plot(xy[0], xy[1], color=c, ls=l)
    for vr, ir, c in zip(vrests, irs, cs):
        ax.plot(vr, ir, 'o', markersize=10, fillstyle='none', color=c)

ax = fig.add_subplot(grid[1, 0])
ax.set_xlabel('V (mV)')
ax.set_ylabel('IK1 (A/F)')
ax.set_xlim(-112, -32)
ax.set_ylim(-0.8, 1.3)
ik1_iv(ax, model2, protocol, ks, vr2)
ax.text(0.65, 0.05, 'After 1 second', transform=ax.transAxes)

ax = fig.add_subplot(grid[1, 1])
ax.set_xlabel('V (mV)')
ax.set_ylabel('IK1 (A/F)')
ax.set_xlim(-112, -32)
ax.set_ylim(-0.4, 0.4)
ik1_iv(ax, model1, protocol, ks, vr1)
ax.text(0.65, 0.05, 'After 1 second', transform=ax.transAxes)
coords = [
    ((-58, 0.02), (-42, -0.15)),
    ((-66, 0.08), (-62, -0.20)),
    ((-75, 0.22), (-68, 0.35)),
    ((-75, 0.22), (-85, 0.35)),
    ((-70, 0.08), (-74, -0.25)),
]
for k, c in zip(ks, coords):
    ax.annotate(
        f'Vr({k})', xy=c[0], xytext=c[1], textcoords='data',
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3'))


#
# C. INaK
#
def inak(ax, model, protocol, levels, colors):

    # Get some variables
    ko = model.labelx('K_o')

    # Create fixed Ko model for IV curves
    fixed = model.clone()
    fko = fixed.labelx('K_o')
    finak = fixed.labelx('I_NaK')
    fixed_states = list(fixed.states())
    for state in fixed_states:
        if state.label() != 'membrane_potential':
            shared.demote(state)

    # Perform simulations, calculate IVs
    vs = np.linspace(-120, 20, 141)
    ivs = []
    for k in levels:
        print(f'Simulating {k}mM with {model.name()}')

        # Create simulation, set [K]o and do first pre-pacing
        s = myokit.Simulation(model, protocol)
        s.set_tolerance(1e-8, 1e-8)
        s.set_constant(ko, k)
        if pre_time:
            s.run(pre_time, log=[])

        # Update fixed model variables and calculate IV
        fko.set_rhs(k)
        x = s.state()
        for i, state in enumerate(fixed_states):
            if state.label() == 'membrane_potential':
                state.set_state_value(x[i])
            else:
                state.set_rhs(x[i])
        ivs.append((vs, finak.pyfunc()(vs)))

    # Plot
    ax.axhline(0, color='#bbbbbb')
    ax.grid(True)
    for xy, c, l in zip(ivs, colors, ls):
        ax.plot(xy[0], xy[1], color=c, ls=l)

ax = fig.add_subplot(grid[2, 0])
ax.set_xlabel('V (mV)')
ax.set_ylabel('INaK (A/F)')
ax.set_xlim(-102, 22)
ax.set_ylim(0.00, 0.36)
inak(ax, model2, protocol, ks, cs)
ax.text(0.65, 0.05, 'After 1 second', transform=ax.transAxes)

ax = fig.add_subplot(grid[2, 1])
ax.set_xlabel('V (mV)')
ax.set_ylabel('INaK (A/F)')
ax.set_xlim(-102, 22)
ax.set_ylim(0.00, 0.36)
inak(ax, model1, protocol, ks, cs)
ax.text(0.65, 0.05, 'After 1 second', transform=ax.transAxes)


# Align labels
fig.align_labels()

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
