#!/usr/bin/env python3
#
# Model figure 2: Instantaneous IV curves at low potassium
#
# Panel A:
#   AP with labels
# Panel B:
#   Quasi-instantaneous IV curves: total current
# Panels C, D, E:
#   Quasi-instantaneous IV curves: IK1, IKur, Ito
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
fname = 'figure-m2-iv'
debug = 'debug' in sys.argv

# Load model
model1 = shared.model()
name = model1.name()

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 800

# Time to pre-pace after changing [K]o level
pre_time = 1000

# Times to evaluate at
if 'courtemanche' in name.lower():
    times = [90, 175, 275]
elif 'grandi' in name.lower():
    times = [100, 225, 350]
elif 'koivumaki' in name.lower():
    times = [100, 225, 350]
elif 'maleckar' in name.lower():
    times = [90, 160, 220]
elif 'ni' in name.lower():
    times = [100, 200, 275]
elif 'nygren' in name.lower():
    times = [100, 210, 320]
elif 'voigt' in name.lower():
    times = [100, 210, 320]
else:
    raise NotImplementedError('No times chosen for this model.')

# Prepare model
shared.prepare_model(model1, protocol, fix_cleft_ko=True, pre_pace=not debug)
model2 = model1.clone()

# Get some variables
t1 = model1.timex()
v1 = model1.labelx('membrane_potential')
ko1 = model1.labelx('K_o')
t2 = model2.timex()
v2 = model2.labelx('membrane_potential')
ko2 = model2.labelx('K_o')

# Set Ko
k1 = 5.4
ko1.set_rhs(k1)
k2 = 2.5
ko2.set_rhs(k2)

# Simulate an AP
s1 = myokit.Simulation(model1, protocol)
s2 = myokit.Simulation(model2, protocol)
s1.set_tolerance(1e-8, 1e-8)
s2.set_tolerance(1e-8, 1e-8)
s1.pre(1000) # 1 second pre-pacing: Update plotting code if changed!
s2.pre(1000) # 1 second pre-pacing: Update plotting code if changed!
d1 = s1.run(tmax).npview()
d2 = s2.run(tmax).npview()


#
# Create figure
#
fig = plt.figure(figsize=(9, 5))  # Two-column size
fig.subplots_adjust(0.08, 0.10, 0.95, 0.98)
grid = GridSpec(2, 3, wspace=0.45, hspace=0.45)

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
fig.text(0.004, 0.95, 'A', letter_font)
fig.text(0.322, 0.95, 'B', letter_font)
fig.text(0.004, 0.43, 'C', letter_font)
fig.text(0.332, 0.43, 'D', letter_font)
fig.text(0.660, 0.43, 'E', letter_font)

# A: AP waveform with labels
ax1 = fig.add_subplot(grid[0, 0])
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('V (mV)')
ax1.plot(pre_time + d1.time(), d1[v1], '-', label=f'{k1} mM', color='k')
ax1.plot(pre_time + d2.time(), d2[v2], '--', label=f'{k2} mM', color='k')
ax1.legend(loc='upper right')
for time in times:
    ax1.plot(pre_time + time, d1.interpolate_at('membrane.V', time), 's')
plt.gca().set_prop_cycle(None)
for time in times:
    ax1.plot(pre_time + time, d2.interpolate_at('membrane.V', time), 'o')


# B: Instantaneous IV for membrane current
ax2 = fig.add_subplot(grid[0, 1:])
ax2.set_xlabel('V (mV)')
ax2.set_ylabel('I (A/F)')
ax2.axhline(0, color='#cccccc')
ax2.set_xlim(-110, 50)
ax2.set_ylim(-1.5, 1.5)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model1, protocol, times, 'cellular_current', s1)
for i, voltage in enumerate(vts):
    ax2.plot(iv_voltages, iv_currents[i])
    ax2.plot(voltage, its[i], 'ks')
plt.gca().set_prop_cycle(None)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model2, protocol, times, 'cellular_current', s2)
for i, voltage in enumerate(vts):
    ax2.plot(iv_voltages, iv_currents[i], '--')
    ax2.plot(voltage, its[i], 'ko')


# C: Instantaneous IV for IK1
ax4 = fig.add_subplot(grid[1, 0])
ax4.set_xlabel('V (mV)')
ax4.set_ylabel('IK1 (A/F)')
ax4.axhline(0, color='#cccccc')
ax4.set_xlim(-110, 50)
ax4.set_ylim(-0.2, 0.9)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model1, protocol, times, 'I_K1', s1)
for i, voltage in enumerate(vts):
    ax4.plot(iv_voltages, iv_currents[i])
    ax4.plot(voltage, its[i], 'ks')
plt.gca().set_prop_cycle(None)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model2, protocol, times, 'I_K1', s2)
for i, voltage in enumerate(vts):
    ax4.plot(iv_voltages, iv_currents[i], '--')
    ax4.plot(voltage, its[i], 'ko')


# D: Instantaneous IV for IKur
ax5 = fig.add_subplot(grid[1, 1])
ax5.set_xlabel('V (mV)')
ax5.set_ylabel('IKur (A/F)')
ax5.axhline(0, color='#cccccc')
ax5.set_xlim(-110, 50)
ax5.set_ylim(-0.15, 0.7)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model1, protocol, times, 'I_Kur', s1)
for i, voltage in enumerate(vts):
    ax5.plot(iv_voltages, iv_currents[i])
    ax5.plot(voltage, its[i], 'ks')
plt.gca().set_prop_cycle(None)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model2, protocol, times, 'I_Kur', s2)
for i, voltage in enumerate(vts):
    ax5.plot(iv_voltages, iv_currents[i], '--')
    ax5.plot(voltage, its[i], 'ko')


# E: Instantaneous IV for IKr
ax6 = fig.add_subplot(grid[1, 2])
ax6.set_xlabel('V (mV)')
ax6.set_ylabel('Ito (A/F)')
ax6.axhline(0, color='#cccccc')
ax6.set_xlim(-100, 50)
ax6.set_ylim(-0.15, 0.7)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model1, protocol, times, 'I_to', s1)
for i, voltage in enumerate(vts):
    ax6.plot(iv_voltages, iv_currents[i])
    ax6.plot(voltage, its[i], 'ks')
plt.gca().set_prop_cycle(None)
vts, its, iv_voltages, iv_currents = shared.quasi_instantaneous_iv(
    model2, protocol, times, 'I_to', s2)
for i, voltage in enumerate(vts):
    ax6.plot(iv_voltages, iv_currents[i], '--')
    ax6.plot(voltage, its[i], 'ko')


# Show / store
path = os.path.join(fpath, f'{fname}-{name}')
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
