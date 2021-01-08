#!/usr/bin/env python3
#
# Membrane resistance
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-a5-rm'

# Load model
model = shared.model('voigt')
name = model.name()

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 100

# Time to pre-pace after changing [K]o level
pre_time = 1000

k = 5.4
dt = 5
dv = 1


# Prepare model
shared.prepare_model(model, protocol, fix_cleft_ko=True, pre_pace=False)

# Get some variables
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')
i_ion_na = model.labelx('cellular_current_nA')
i_ion = model.labelx('cellular_current')
clamp = model.labelx('clamp')

# Simulate normal AP
print(f'Simulating {k}mM')
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
s.set_constant(ko, k)
if pre_time:
    s.pre(pre_time)
d = s.run(tmax).npview().npview()

# Create figure
fig = plt.figure(figsize=(9, 10))  # Two-column size
fig.subplots_adjust(0.11, 0.07, 0.985, 0.99)
#grid = matplotlib.gridspec.GridSpec(2, nk, wspace=0.04, hspace=0.04)

ax1 = fig.add_subplot(4, 1, 1)
ax2 = fig.add_subplot(4, 1, 2)
ax3 = fig.add_subplot(2, 1, 2)

ax1.set_ylabel('V (mV)')
ax2.set_ylabel('I net (A/F)')
ax3.set_ylabel('I (A/F)')
ax3.set_xlabel('Time (ms)')

# Plot normal AP
ax1.plot(1e3 + d.time(), d[vm], ':', color='#999999')
ax2.plot(1e3 + d.time(), d[i_ion], ':', color='#999999')

color0 = 'black'
color1 = 'black'
color2 = 'red'

# Run until time t
t = 25
s.reset()
s.set_constant(clamp, 0)
d0 = s.run(t).npview()
x = s.state()
ax1.plot(1e3 + d0.time(), d0[vm], color=color0)
ax2.plot(1e3 + d0.time(), d0[i_ion], color=color0)

#currents = shared.currents(model)
#currents.remove('ipca.IpCa')
#currents.remove('ikur.IKur')
#currents.remove('ikr.IKr')
#currents.remove('iclca.IClCa')
#currents.remove('iks.IKs')
#currents.remove('ikp.IKp')
#currents.remove('ikach.IKACh')
#currents.remove('inal.INaL')
#currents.sort(key=lambda x: -d0[x][0])
currents = [
    'ik1.IK1',
    'inak.INaK',
    'ito.Ito',
    'ina.INa',
    'ical.ICaL',
    'inaca.INaCa',
    'inab.INaB',
    'icab.ICaB',
    'iclb.IClB',
]
cmap = matplotlib.cm.get_cmap('tab10')
for i, current in enumerate(currents):
    label = current[current.index('.') + 1:]
    ax3.plot(1e3 + d0.time(), d0[current], label=label, color=cmap(i))

# Run with +dv for dt time units
x1 = list(x)
x1[vm.indice()] += dv
s.set_state(x1)
s.set_constant(clamp, 1)

d1 = s.run(dt).npview()
i1 = d1[i_ion_na]
ax1.plot(1e3 + d1.time(), d1[vm], color=color1)
ax2.plot(1e3 + d1.time(), d1[i_ion], color=color1)
for i, current in enumerate(currents):
    ax3.plot(1e3 + d1.time(), d1[current], color=cmap(i))

p = {'alpha': 0.2, 'lw': 2, 'ls': '-'}
ax1.plot([1e3 + t, 1e3 + t], [d0[vm][-1], d1[vm][0]], color=color1, **p)
ax2.plot([1e3 + t, 1e3 + t], [d0[i_ion][-1], d1[i_ion][0]], color=color1, **p)
for i, c in enumerate(currents):
    ax3.plot([1e3 + t, 1e3 + t], [d0[c][-1], d1[c][0]], color=cmap(i), **p)

s.set_constant(clamp, 0)
if dt > 1e-3:
    e1 = s.run(49 - s.time()).npview()
    ax1.plot(1e3 + e1.time(), e1[vm], '-', color=color1, alpha=0.2)
    ax2.plot(1e3 + e1.time(), e1[i_ion], '-', color=color1, alpha=0.2)
    for i, current in enumerate(currents):
        ax3.plot(1e3 + e1.time(), e1[current], color=cmap(i), alpha=0.2)


# Run with -dv for dt time units
x2 = list(x)
x2[vm.indice()] -= dv
s.set_state(x2)
s.set_time(t)
s.set_constant(clamp, 1)

d2 = s.run(dt).npview()
i2 = d2[i_ion_na]
ax1.plot(1e3 + d2.time(), d2[vm], '--', color=color2)
ax2.plot(1e3 + d2.time(), d2[i_ion], '--', color=color2)
for i, current in enumerate(currents):
    ax3.plot(1e3 + d2.time(), d2[current], '--', color=cmap(i))

p = {'alpha': 0.2, 'lw': 2, 'ls': '-'}
ax1.plot([1e3 + t, 1e3 + t], [d0[vm][-1], d2[vm][0]], color=color2, **p)
ax2.plot([1e3 + t, 1e3 + t], [d0[i_ion][-1], d2[i_ion][0]], color=color2, **p)
for i, c in enumerate(currents):
    ax3.plot([1e3 + t, 1e3 + t], [d0[c][-1], d2[c][0]], color=cmap(i), **p)

s.set_constant(clamp, 0)
if dt > 1e-3:
    e2 = s.run(49 - s.time()).npview()
    ax1.plot(1e3 + e2.time(), e2[vm], '--', color=color2, alpha=0.2)
    ax2.plot(1e3 + e2.time(), e2[i_ion], '--', color=color2, alpha=0.2)
    for i, current in enumerate(currents):
        ax3.plot(1e3 + e2.time(), e2[current], '--', color=cmap(i), alpha=0.2)

# Set padding on fig 2
lo = min(np.min(d1[i_ion]), np.min(d2[i_ion]))
hi = max(np.max(d1[i_ion]), np.max(d2[i_ion]))
pad = (hi - lo) * 0.04
ax2.set_ylim(lo - pad, hi + pad)

# Calculate R
r = 2 * dv / (i1[-1] - i2[-1])
ax2.text(0.8, 0.6, f'Rm = {round(r, 1)}', transform=ax2.transAxes)

if dt > 1e-3:
    xlim = 999, 1101
else:
    xlim = 1025 - 2 * dt, 1025 + 2 * dt

ax1.set_xlim(*xlim)
ax2.set_xlim(*xlim)
ax3.set_xlim(*xlim)

ax3.legend(loc='center right')

ax1.text(0.02, 0.9, f'dt = {dt}', transform=ax1.transAxes)
ax1.text(0.02, 0.8, f'dv = {dv}', transform=ax1.transAxes)
ax1.text(0.02, 0.7, f'[K]o = {k}', transform=ax1.transAxes)

if dt > 1e-3:
    props = {'arrowprops': {'arrowstyle': '->'}, 'horizontalalignment': 'center'}
    y = np.max(d0[vm])
    ax1.annotate('Unclamped', (1010, y + 1), (1010, y + 15), **props)
    y = np.max(d1[vm])
    ax1.annotate('Clamped', (1027, y + 1), (1027, y + 15), **props)
    y = e1.interpolate_at(vm, 40)
    ax1.annotate('Unclamped', (1040, y + 1), (1040, y + 15), **props)

fig.align_labels()

# Show / store
path = f'{fname}-ko-{k}-dv-{dv}-dt-{dt}'
path = os.path.join(fpath, path.replace('.', '_'))
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
