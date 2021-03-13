#!/usr/bin/env python3
#
# Short and long-term changes in INaK
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
fname = 'figure-s7-fast-slow'
debug = 'debug' in sys.argv

# Load model
model = shared.model()
name = model.name()
shared.splash(fname, name)

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 800

# Time to pre-pace after changing [K]o level
pre_time_1 = 1000
pre_time_2 = 15 * 60 * 1000

# Prepare model
shared.prepare_model(model, protocol, pre_pace=not debug)

# Get some variables
v = model.labelx('membrane_potential')
t = model.time()
ko = model.labelx('K_o')
ki = model.labelx('K_i')
nai = model.labelx('Na_i')
cai = model.label('Ca_i')

# Create fixed Ko model for IV curves
fixed = model.clone()
fki = fixed.labelx('K_i')
fko = fixed.labelx('K_o')
fik1 = fixed.labelx('I_K1')
finak = fixed.labelx('I_NaK')
finaca = fixed.labelx('I_NaCa')
fixed_states = list(fixed.states())
for state in fixed_states:
    if state.label() != 'membrane_potential':
        shared.demote(state)


# Perform simulations, calculate IVs
ks = shared.ko_levels
cs = shared.ko_colors
d1s = []
d2s = []
d3s = []
vs = np.linspace(-120, 20, 141)
ik1_1, ik1_2 = [], []
inak_1, inak_2 = [], []
inaca_1, inaca_2 = [], []
for k in ks:
    print(f'Simulating {k}mM')

    # Create simulation, set [K]o and do first pre-pacing
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    s.set_constant(ko, k)
    if pre_time_1:
        s.run(pre_time_1, log=[])

    # Update model variables and calculate early IV
    fko.set_rhs(k)
    x = s.state()
    for i, state in enumerate(fixed_states):
        if state.label() == 'membrane_potential':
            state.set_state_value(x[i])
        else:
            state.set_rhs(x[i])
    ik1_1.append((vs, fik1.pyfunc()(vs)))
    inak_1.append((vs, finak.pyfunc()(vs)))
    inaca_1.append((vs, finaca.pyfunc()(vs)))

    # Simulate early AP
    d1s.append(s.run(tmax).npview())
    if debug:
        continue

    # Second pre-pacing
    s.reset()
    log_vars = [t, ki, nai]
    if cai is not None:
        log_vars.append(cai)
    d3s.append(s.run(
        pre_time_2, log=log_vars, log_interval=cl).npview())

    # Update model variables and calculate late IV
    x = s.state()
    for i, state in enumerate(fixed_states):
        if state.label() == 'membrane_potential':
            state.set_state_value(x[i])
        else:
            state.set_rhs(x[i])
    ik1_2.append((vs, fik1.pyfunc()(vs)))
    inak_2.append((vs, finak.pyfunc()(vs)))
    inaca_2.append((vs, finaca.pyfunc()(vs)))

    # Simulate late AP
    d2s.append(s.run(tmax).npview())
    del(s)
if debug:
    d3s = d2s = d1s
    ik1_2 = ik1_1
    inak_2 = inak_1
    inaca_2 = inaca_1


#
# Create figure
#
fig = plt.figure(figsize=(9, 11))  # Two-column size
fig.subplots_adjust(0.09, 0.05, 0.95, 0.98)
grid = matplotlib.gridspec.GridSpec(5, 6, wspace=1.40, hspace=0.5)

# Add model name
fig.patches.append(matplotlib.patches.Rectangle(
    (0.965, 0), 0.035, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.975, 0.5, name, name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.001, 0.96, 'A', letter_font)
fig.text(0.322, 0.96, 'B', letter_font)
if cai is not None:
    fig.text(0.655, 0.96, 'C', letter_font)
fig.text(0.001, 0.78, 'D', letter_font)
fig.text(0.485, 0.77, 'E', letter_font)
fig.text(0.001, 0.56, 'F', letter_font)
fig.text(0.485, 0.56, 'G', letter_font)
#fig.text(0.001, 0.37, 'H', letter_font)
#fig.text(0.485, 0.37, 'I', letter_font)
#fig.text(0.001, 0.17, 'J', letter_font)
#fig.text(0.485, 0.17, 'K', letter_font)

# A: Ki over time
ax = fig.add_subplot(grid[0, :2])
ax.set_xlabel('Time (s)')
ax.set_ylabel('[K]i (mM)')
for d, k, c in zip(d3s, ks, cs):
    ax.plot(d.time() * 1e-3, d[ki], color=c)

# B: Nai over time
ax = fig.add_subplot(grid[0, 2:4])
ax.set_xlabel('Time (s)')
ax.set_ylabel('[Na]i (mM)')
for d, k, c in zip(d3s, ks, cs):
    ax.plot(d.time() * 1e-3, d[nai], color=c)

# C: Cai over time
if cai is not None:
    ax = fig.add_subplot(grid[0, 4:])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('[Ca]i (uM)')
    for d, k, c in zip(d3s, ks, cs):
        ax.plot(d.time() * 1e-3, d[cai] * 1e3, color=c)

# D: Immediate change in AP waveform
ax = fig.add_subplot(grid[1, :3])
ax.set_xlabel('Time (s)')
ax.set_ylabel('V (mV)')
for d, k, c in zip(d1s, ks, cs):
    ax.plot(d.time() * 1e-3, d[v], label=f'{k}', color=c)
ax.legend(loc='upper right')

# E: Late change in AP waveform
ax = fig.add_subplot(grid[1, 3:])
ax.set_xlabel('Time (s)')
ax.set_ylabel('V (mV)')
for d, k, c in zip(d2s, ks, cs):
    ax.plot(d.time() * 1e-3, d[v], color=c)

# F: IK1 IV early
ax = fig.add_subplot(grid[2:, :3])
ax.set_xlabel('V (mV)')
ax.set_ylabel('I (A/F)')
#ax.grid(True)
ax.axhline(0, color='#bbbbbb')
label = True
for x, y, z, c in zip(ik1_1, inak_1, inaca_1, cs):
    ax.plot(x[0], x[1], color=c, label='IK1' if label else None)
    ax.plot(y[0], y[1], '--', color=c, label='INaK' if label else None)
    ax.plot(z[0], z[1], ':', color=c, label='INaCa' if label else None)
    label = False
ax.set_ylim(-0.5, 0.90)
ax.text(0.043, 0.96, 'After 1 second', transform=ax.transAxes)
ax.legend(loc='lower right')

# G: IK1 IV late
ax = fig.add_subplot(grid[2:, 3:])
ax.set_xlabel('V (mV)')
ax.set_ylabel('I (A/F)')
ax.axhline(0, color='#bbbbbb')
#ax.grid(True)
for x, y, z, c in zip(ik1_2, inak_2, inaca_2, cs):
    ax.plot(x[0], x[1], color=c)
    ax.plot(y[0], y[1], '--', color=c)
    ax.plot(z[0], z[1], ':', color=c)
ax.set_ylim(-0.5, 0.90)
ax.text(0.043, 0.96, 'After 15 minutes', transform=ax.transAxes)


fig.align_labels()

# Show / store
path = os.path.join(fpath, fname + '-' + name)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')

