#!/usr/bin/env python3
#
# Model figure 4: INa availability
#
# Panel A: After 1 second
# Panel B: After 15 minutes
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-m4-ina'
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
pre_time_1 = 1000
pre_time_2 = 15 * 60 * 1000

# Prepare model
shared.prepare_model(model, protocol, fix_cleft_ko=True, pre_pace=not debug)

# Get some variables
t = model.time()
v = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Get steady-state inactivation curve
try:
    # Courtemanche and grandi families use (h * j) for inactivation
    ina_j_inf = model.labelx('ina_j_inf')
    ina_h_inf = model.labelx('ina_h_inf')
    fj = ina_j_inf.pyfunc()
    fh = ina_h_inf.pyfunc()
    f = lambda v: fj(v) * fh(v)
except myokit.IncompatibleModelError as e:
    # Nygren family uses (0.9 * h1 + 0.1 * h2) for inactivation
    ina_h12_inf = model.labelx('ina_h12_inf')
    f = ina_h12_inf.pyfunc()

# Calculate inactivation curve
vs = np.linspace(-160, 40, 200)
ss = f(vs)

# Perform simulations
ks = shared.ko_levels
cs = shared.ko_colors
v1s = []
v2s = []
for k in ks:
    print(f'Simulating {k}mM')

    # Create simulation, set [K]o and do first pre-pacing
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    s.set_constant(ko, k)
    if pre_time_1:
        s.run(pre_time_1, log=myokit.LOG_NONE)

    # Get resting potential
    v1s.append(s.state()[v.indice()])

    # Stop if in a hurry
    if debug:
        v1s = v1s * len(ks)
        v2s = v1s
        break

    # Second pre-pacing
    s.reset()
    s.run(pre_time_2 - pre_time_1, log=myokit.LOG_NONE)

    # Get resting potential
    v2s.append(s.state()[v.indice()])


#
# Create figure
#
fig = plt.figure(figsize=(4.3, 2.0))  # One-column size
fig.subplots_adjust(0.13, 0.22, 0.91, 0.98)
grid = GridSpec(1, 2, wspace=0.05)

# Add model name
fig.patches.append(Rectangle((0.93, 0), 0.07, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.95, 0.5, name, name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
#fig.text(0.004, 0.95, 'A', letter_font)
#fig.text(0.245, 0.95, 'B', letter_font)
#fig.text(0.485, 0.95, 'C', letter_font)

# A: Change after 1 second
ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel('Voltage (mV)')
#ax.set_ylabel('Recovered INa')
ax.set_xlim(-130, -20)
ax.set_ylim(-0.05, 1.25)
ax.plot(vs, ss, color='#999999')
for vr, k, c in zip(v1s, ks, cs):
    ax.plot([vr], [f(vr)], 'o', markersize=8, fillstyle='none', color=c,
            label=f'{k} mM')
ax.text(0.043, 0.88, 'After 1 second', transform=ax.transAxes)

# B: Change after 15 minutes
ax = fig.add_subplot(grid[0, 1])
ax.set_xlabel('Voltage (mV)')
ax.set_ylabel('Recovered INa')
ax.set_yticklabels([])
ax.set_xlim(-130, -20)
ax.set_ylim(-0.05, 1.25)
ax.plot(vs, ss, color='#999999')
for vr, k, c in zip(v2s, ks, cs):
    ax.plot([vr], [f(vr)], 'o', markersize=8, fillstyle='none', color=c,
            label=f'{k}')
ax.text(0.043, 0.88, 'After 15 minutes', transform=ax.transAxes)
ax.legend(loc='center right', handlelength=1, handletextpad=0.5)

# Show / store
path = os.path.join(fpath, fname + '-' + name)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')

