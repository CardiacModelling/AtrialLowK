#!/usr/bin/env python3
#
# Figure 5: Short and long-term changes to INaK and IK1
#
import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import shared
import sys

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-4'
debug = 'debug' in sys.argv

# Load model
model = shared.model('voigt')
name = model.name()

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time and voltage range to show in plots
tmax = 800
vlim = -102, 2

# Time to pre-pace after changing [K]o level
pre_time_1 = 1000
pre_time_2 = 15 * 60 * 1000

# Prepare model
shared.prepare_model(model, protocol, fix_cleft_ko=True, pre_pace=not debug)

# Get some variables
v = model.labelx('membrane_potential')
t = model.time()
ko = model.labelx('K_o')
ki = model.labelx('K_i')
nai = model.labelx('Na_i')

# Create fixed Ko model for IV curves
fixed = model.clone()
fki = fixed.labelx('K_i')
fko = fixed.labelx('K_o')
finak = fixed.labelx('I_NaK')
fik1 = fixed.labelx('I_K1')
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
inak_v1 = []
inak_v2 = []
ik1_v1 = []
ik1_v2 = []
vrs1 = []
vrs2 = []
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
for k in ks:
    print(f'Simulating {k}mM')

    # Create simulation, set [K]o and do first pre-pacing
    s.reset()
    s.set_constant(ko, k)
    if pre_time_1:
        s.run(pre_time_1, log=[])

    # Update model variables and calculate early IVs
    fko.set_rhs(k)
    x = s.state()
    for i, state in enumerate(fixed_states):
        if state.label() == 'membrane_potential':
            state.set_state_value(x[i])
        else:
            state.set_rhs(x[i])
    inak_v1.append((vs, finak.pyfunc()(vs)))
    f = fik1.pyfunc()
    ik1_v1.append((vs, f(vs)))

    # Store vrest and ik1(vrest)
    vr = s.state()[v.indice()]
    vrs1.append((vr, f(vr)))

    # Simulate early AP
    d1s.append(s.run(tmax).npview())
    if debug:
        continue

    # Second pre-pacing
    s.reset()
    d3s.append(s.run(pre_time_2, log=[t, ki, nai]).npview())

    # Update model variables and calculate late IVs
    x = s.state()
    for i, state in enumerate(fixed_states):
        if state.label() == 'membrane_potential':
            state.set_state_value(x[i])
        else:
            state.set_rhs(x[i])
    inak_v2.append((vs, finak.pyfunc()(vs)))
    f = fik1.pyfunc()
    ik1_v2.append((vs, f(vs)))

    # Store vrest and ik1(vrest)
    vr = s.state()[v.indice()]
    vrs2.append((vr, f(vr)))

    # Simulate late AP
    d2s.append(s.run(tmax).npview())
if debug:
    d3s = d2s = d1s
    inak_v2 = inak_v1
    ik1_v2 = ik1_v1
    vrs2 = vrs1

#
# Create figure
#
fig = plt.figure(figsize=(9, 5.5))  # Two-column size
fig.subplots_adjust(0.075, 0.085, 0.99, 0.88)
grid = matplotlib.gridspec.GridSpec(5, 4, wspace=0.60, hspace=1.2)

# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.975, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.002, 0.85, 'A', letter_font)
fig.text(0.505, 0.85, 'B', letter_font)
fig.text(0.002, 0.49, 'C', letter_font)
fig.text(0.505, 0.49, 'D', letter_font)
#fig.text(0.255, 0.91, 'B', letter_font)
#fig.text(0.505, 0.91, 'C', letter_font)
#fig.text(0.755, 0.91, 'D', letter_font)
#fig.text(0.002, 0.52, 'E', letter_font)
#fig.text(0.505, 0.52, 'F', letter_font)

# A: Immediate change in AP waveform
ax = fig.add_subplot(grid[0:2, 0:2])
ax.set_xlabel('Time (s)')
ax.set_ylabel('V (mV)')
for d, k, c in zip(d1s, ks, cs):
    ax.plot(d.time() * 1e-3, d[v], label=f'{k} mM', color=c)
ax.legend(loc=(0.352, 1.07), ncol=len(ks))
ax.text(0.710, 0.87, 'After 1 second', transform=ax.transAxes)

# D: Late change in AP waveform
ax = fig.add_subplot(grid[0:2, 2:4])
ax.set_xlabel('Time (s)')
ax.set_ylabel('V (mV)')
for d, k, c in zip(d2s, ks, cs):
    ax.plot(d.time() * 1e-3, d[v], color=c)
ax.text(0.665, 0.87, 'After 15 minutes', transform=ax.transAxes)

# E: IK1 and INaK IV early
ax = fig.add_subplot(grid[2:, :2])
ax.set_xlabel('V (mV)')
ax.set_ylabel('I (A/F)')
ax.set_xlim(*vlim)
ax.set_ylim(-0.12, 0.39)
ax.axhline(0, color='#bbbbbb')


for k, xy, qr, c in zip(ks, ik1_v1, vrs1, cs):
    label = 'IK1' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, label=label)
    ax.plot(qr[0], qr[1], 'o', color=c, fillstyle='none')
    ax.axvline(qr[0], color=c, alpha=0.75, lw=1, zorder=0)
ax.annotate(
    'Vrest',
    xy=(vrs1[0][0], vrs1[0][1] - 0.01),
    xytext=(vrs1[0][0] + 3, -0.04),
    textcoords='data',
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3'))
for k, xy, c in zip(ks, inak_v1, cs):
    label = 'INaK' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, ls='--', label=label)
ax.text(0.710, 0.92, 'After 1 second', transform=ax.transAxes)
ax.legend(ncol=2, loc='lower right')
#ax.grid(True)

# F: IK1 and INaK IV late
ax = fig.add_subplot(grid[2:, 2:])
ax.set_xlabel('V (mV)')
ax.set_ylabel('IK1 (A/F)')
ax.set_xlim(*vlim)
ax.set_ylim(-0.12, 0.39)
ax.axhline(0, color='#bbbbbb')
for k, xy, qr, c in zip(ks, ik1_v2, vrs2, cs):
    label = 'IK1' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, label=label)
    ax.plot(qr[0], qr[1], 'o', color=c, fillstyle='none')
    ax.axvline(qr[0], color=c, alpha=0.75, lw=1, zorder=0)
ax.text(0.665, 0.92, 'After 15 minutes', transform=ax.transAxes)
for k, xy, c in zip(ks, inak_v2, cs):
    label = 'INaK' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, ls='--', label=label)
ax.legend(ncol=2, loc='lower right')
#ax.grid(True)

fig.align_labels()

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')

