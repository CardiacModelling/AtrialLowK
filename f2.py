#!/usr/bin/env python3
#
# Figure 2: Currents in Grandi model and modified model.
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import numpy as np
import os
import shared
import sys

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-2'
debug = 'debug' in sys.argv

# Load model
model1 = shared.model('grandi')
model2 = shared.model('voigt')

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 600 + 50

# Prepare models
shared.prepare_model(model1, protocol, pre_pace=not debug)
shared.prepare_model(model2, protocol, pre_pace=not debug)

# Variables to plot, and y limits
conc = {
    'potassium.K_i': '$[K^+]_i$ (mM)',
    'sodium.Na_i': '$[Na^+]_i$ (mM)',
    'chloride.Cl_i': '$[Cl^-]_i$ (mM)',
}
variables = {
    'ina.INa': 'INa (A/F)',
    'ical.ICaL': 'ICaL (A/F)',
    'ikur.IKur': 'IKur (A/F)',
    'ito.Ito': 'Ito (A/F)',
    'ikr.IKr': 'IKr (A/F)',
    'iks.IKs': 'IKs (A/F)',
    'calcium.Ca_i': '$[Ca^{2+}]_i$ (μM)',
    'inab.INaB': 'INaB (A/F)',
    'icab.ICaB': 'ICaB (A/F)',
    'iclb.IClB': 'IClB (A/F)',
    'ik1.IK1': 'IK1 (A/F)',
    'inak.INaK': 'INaK (A/F)',
    'inaca.INaCa': 'INaCa (A/F)',
    'ikach.IKACh': 'IKACh (A/F)',

}
changed = set(['ina.INa', 'ik1.IK1', 'ikach.IKACh'])


# Run simulation
print(f'Simulating...')
s1 = myokit.Simulation(model1, protocol)
s1.set_tolerance(1e-8, 1e-8)
d1 = s1.run(tmax).npview()
s2 = myokit.Simulation(model2, protocol)
s2.set_tolerance(1e-8, 1e-8)
d2 = s2.run(tmax).npview()

# Show without offset
d1['engine.time'] -= 50
d2['engine.time'] -= 50

# Times to show
xlim = (-20, 520)

#
# Create figure
#
fig = plt.figure(figsize=(9, 9))  # Two-column figure
fig.subplots_adjust(0.12, 0.05, 0.995, 0.99)
n = (1 + len(variables)) // 2
grid = matplotlib.gridspec.GridSpec(2 + n, 2, hspace=0.085, wspace=0.4)

# Add panel letters
letter_font = {
    'weight': 'bold', 'fontsize': 16, 'horizontalalignment': 'center'}
fig.text(0.013, 0.980, 'A', letter_font)
fig.text(0.013, 0.760, 'B', letter_font)
fig.text(0.013, 0.655, 'C', letter_font)
fig.text(0.013, 0.550, 'D', letter_font)
fig.text(0.013, 0.445, 'E', letter_font)
fig.text(0.013, 0.340, 'F', letter_font)
fig.text(0.013, 0.235, 'G', letter_font)
fig.text(0.013, 0.130, 'H', letter_font)

fig.text(0.520, 0.980, "A'", letter_font)
fig.text(0.520, 0.760, 'I', letter_font)
fig.text(0.520, 0.655, 'J', letter_font)
fig.text(0.520, 0.550, 'K', letter_font)
fig.text(0.520, 0.445, 'L', letter_font)
fig.text(0.520, 0.340, 'M', letter_font)
fig.text(0.520, 0.235, 'N', letter_font)
fig.text(0.520, 0.130, 'O', letter_font)

# Action potential
ax = fig.add_subplot(grid[0:2, 0])
ax.set_xlabel('Time (ms)')
ax.set_ylabel('V (mV)')
ax.set_xlim(*xlim)
ax.plot(d1.time(), d1['membrane.V'], label=shared.fancy_name(model1))
ax.plot(d2.time(), d2['membrane.V'], label=shared.fancy_name(model2))
ax.set_xticklabels([])
ax.set_yticks([-80, -40, 0, 40])
ax.legend()

vr1 = np.round(d1['membrane.V'][0], 1)
vr2 = np.round(d2['membrane.V'][0], 1)

ax.text(0.7, 0.48, f'Vr = {vr1} mV',  transform=ax.transAxes,
        horizontalalignment='left', color='tab:blue')
ax.text(0.7, 0.37, f'Vr = {vr2} mV',  transform=ax.transAxes,
        horizontalalignment='left', color='tab:orange')

ax = fig.add_subplot(grid[0:2, 1])
ax.set_xlabel('Time (ms)')
ax.set_ylabel('V (mV)')
ax.set_xlim(*xlim)
ax.plot(d1.time(), d1['membrane.V'], label=shared.fancy_name(model1))
ax.plot(d2.time(), d2['membrane.V'], label=shared.fancy_name(model2))
ax.set_xticklabels([])
ax.set_yticks([-80, -40, 0, 40])
ax.legend()

# Variables
for i, variable in enumerate(variables):
    label = variables[variable]
    if label is None:
        label = variable.split('.')[1]

    row = 2 + i % n
    col = i // n

    ax = fig.add_subplot(grid[row, col])
    ax.set_xlim(*xlim)
    label_prop = {}
    #factor = factors.get(variable, 1)
    if variable in changed:
        label_prop['fontweight'] = 'bold'
        #if factor != 1:
        #    ax.text(0.9, 0.5, f'x {factor}',  transform=ax.transAxes,
        #            horizontalalignment='right', fontweight='bold')

    ax.set_ylabel(label, **label_prop)
    ax.tick_params('y', labelsize='small')
    if row == n + 1 or (row == n and col == 1):
        ax.set_xlabel('Time (ms)')
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis='x', direction='in')
    try:
        ax.plot(d1.time(), d1[variable])
    except KeyError:
        ax.plot([0], [0])   # For color cycle
    ax.plot(d2.time(), d2[variable])

    # INa inset
    if variable == 'ina.INa':
        axi = ax.inset_axes([0.20, 0.35, 0.5, 0.5])
        axi.tick_params(labelsize='small')
        axi.set_xlim(-0.5, 3)
        axi.plot(d1.time(), d1['ina.INa'])
        axi.plot(d2.time(), d2['ina.INa'])

fig.align_labels()


# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
