#!/usr/bin/env python3
#
# Figure 3: Short and long-term changes to INaK and IK1
#
import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-3'
shared.splash(fname)

# Load model
model = shared.model('voigt')
name = model.name()

# Create protocol
protocol = shared.default_protocol()

# Maximum time and voltage range to show in plots
tmax = 800
vlim = -95, -20

# Time to pre-pace after changing [K]o level
pre_time_1 = 1000
pre_time_2 = 15 * 60 * 1000

# Prepare model
shared.prepare_model(model, protocol)

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
    print(f'Running for {k}mM')

    # Get state after 1 beat
    x = shared.state_fast(model, protocol, k)

    # Update model variables and calculate early IVs
    fko.set_rhs(k)
    for i, state in enumerate(fixed_states):
        if state.label() == 'membrane_potential':
            state.set_state_value(x[i])
        else:
            state.set_rhs(x[i])

    vr = [x[v.indice()]]

    f = finak.pyfunc()
    vr.append(f(vr[0]))
    inak_v1.append((vs, f(vs)))

    f = fik1.pyfunc()
    vr.append(f(vr[0]))
    ik1_v1.append((vs, f(vs)))

    vrs1.append(vr)

    # Get early AP
    d1s.append(shared.log_fast(model, protocol, k))

    # Update to state after 900 beats
    d3s.append(shared.log_progress(model, protocol, k))

    # Update model variables and calculate late IVs
    x = shared.state_slow(model, protocol, k)
    for i, state in enumerate(fixed_states):
        if state.label() == 'membrane_potential':
            state.set_state_value(x[i])
        else:
            state.set_rhs(x[i])

    vr = [x[v.indice()]]

    f = finak.pyfunc()
    vr.append(f(vr[0]))
    inak_v2.append((vs, f(vs)))

    f = fik1.pyfunc()
    vr.append(f(vr[0]))
    ik1_v2.append((vs, f(vs)))

    vrs2.append(vr)

    # Simulate late AP
    d2s.append(shared.log_slow(model, protocol, k))


#
# Create figure
#
fig = plt.figure(figsize=(9, 6.5))  # Two-column size
fig.subplots_adjust(0.085, 0.075, 0.98, 0.95)
grid = matplotlib.gridspec.GridSpec(5, 4, wspace=0.60, hspace=1.0)

# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.975, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.002, 0.930, 'A', letter_font)
fig.text(0.495, 0.920, 'B', letter_font)
fig.text(0.002, 0.546, 'C', letter_font)
fig.text(0.495, 0.546, 'D', letter_font)

# A: Immediate change in AP waveform
ax = fig.add_subplot(grid[0:2, 0:2])
ax.set_xlabel('Time (s)')
ax.set_ylabel('V (mV)')
for d, k, c in zip(d1s, ks, cs):
    ax.plot(d.time() * 1e-3, d[v], label=f'{k} mM', color=c)
#ax.legend(loc=(0.352, 1.07), ncol=len(ks))
ax.legend(loc='upper right', ncol=1)
ax.text(0.15, 0.87, 'After 1 second', transform=ax.transAxes)

# D: Late change in AP waveform
ax = fig.add_subplot(grid[0:2, 2:4])
ax.set_xlabel('Time (s)')
ax.set_ylabel('V (mV)')
for d, k, c in zip(d2s, ks, cs):
    ax.plot(d.time() * 1e-3, d[v], color=c)
ax.text(0.15, 0.87, 'After 15 minutes', transform=ax.transAxes)

# C: IK1 and INaK IV early
ax = fig.add_subplot(grid[2:, :2])
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.set_xlabel('V (mV)')
ax.set_ylabel('I (A/F)')
ax.set_xlim(*vlim)
ax.set_ylim(-0.10, 0.30)
ax.axhline(0, color='#bbbbbb')
for k, xy, qr, c in zip(ks, ik1_v1, vrs1, cs):
    label = '$I_{K1}$' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, label=label)
    ax.plot(qr[0], qr[1], 'o', color=c, fillstyle='none')
    ax.plot(qr[0], qr[2], 'o', color=c, fillstyle='none')
    ax.axvline(qr[0], color=c, lw=0.7, alpha=0.8, zorder=1)
ax.annotate(
    '$V_r$',
    xy=(vrs1[0][0], vrs1[0][2] - 0.01),
    xytext=(vrs1[0][0] + 3, -0.04),
    textcoords='data',
    arrowprops=dict(arrowstyle='->', connectionstyle='arc3'))
for k, xy, c in zip(ks, inak_v1, cs):
    label = '$I_{NaK}$' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, ls='--', label=label)
ax.text(0.05, 0.95, 'After 1 second', transform=ax.transAxes,
        bbox=dict(boxstyle='square', ec='#ffffff', fc='#ffffff'))
ax.legend(ncol=2, loc='lower right')
#ax.grid(True)

# D: IK1 INaK IV late
ax = fig.add_subplot(grid[2:, 2:])
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.set_xlabel('V (mV)')
ax.set_ylabel('I (A/F)')
ax.set_xlim(*vlim)
ax.set_ylim(-0.10, 0.30)
ax.axhline(0, color='#bbbbbb')
for k, xy, qr, c in zip(ks, ik1_v2, vrs2, cs):
    label = '$I_{K1}$' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, label=label)
    ax.plot(qr[0], qr[1], 'o', color=c, fillstyle='none')
    ax.plot(qr[0], qr[2], 'o', color=c, fillstyle='none')
    ax.axvline(qr[0], color=c, lw=0.7, alpha=0.8, zorder=1)
ax.text(0.05, 0.95, 'After 15 minutes', transform=ax.transAxes,
        bbox=dict(boxstyle='square', ec='#ffffff', fc='#ffffff'))
for k, xy, c in zip(ks, inak_v2, cs):
    label = '$I_{NaK}$' if k == 5.4 else None
    ax.plot(xy[0], xy[1], color=c, ls='--', label=label)
ax.legend(ncol=2, loc='lower right')
#ax.grid(True)

fig.align_labels()

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
plt.savefig(path + '.jpg', dpi=600)
print('Done')

