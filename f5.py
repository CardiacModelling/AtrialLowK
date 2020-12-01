#!/usr/bin/env python3
#
# Combine strand simulation results for multiple concentrations
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import numpy as np
import os
import shared
import sys

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-5'
debug = 'debug' in sys.argv

# Settings
af = False
cvt = 40
n_beats = 1
i_beat = 51
i_final = i_beat + 1

# Long or short-term?
if 'long' in sys.argv:
    fname = 'figure-6'
    n_beats = 15 * 60
    i_beat += n_beats - 1
    i_final = i_beat + 1

# Load model
model = shared.model('voigt')
name = model.name()
print(f'Loaded {name}')

# Get data file names
dpath = 'strand-data'
spath = 'sd-data'
name1 = name + ('-af' if af else '')
name2 = f'{name1}-cvt-{cvt}'

# Set af mode
if af:
    print('Enabling AF mode')
    shared.enable_af(model)

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Prepare model
shared.prepare_model(model, protocol, fix_cleft_ko=True, pre_pace=not debug)

# Get some variables
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Ko levels and colours
ks = shared.ko_levels
cs = shared.ko_colors
nk = len(ks)

# Get steady-state inactivation curve
ina_j_inf = model.labelx('ina_j_inf')
ina_h_inf = model.labelx('ina_h_inf')
fj = ina_j_inf.pyfunc()
fh = ina_h_inf.pyfunc()
f = lambda v: fj(v) * fh(v)
ftm = model.get('ina.m.tau').pyfunc()

# Calculate inactivation curve
vs = np.linspace(-160, 40, 200)
ss = f(vs)

# Get Vrest for inf(h,j), and get _actual_ availability
vrs = []
avs = []
ms = []
for k in ks:
    # Load state at start of final low-K pre-pacing beat
    dname = f'{name2}-ko-{str(float(k)).replace(".", "-")}-state{i_beat}.zip'
    path = os.path.join(dpath, dname)
    if os.path.isfile(path):
        state_final = myokit.load_state_bin(path)
        # Get state in central cell
        n = model.count_states()
        state_final = state_final[0*n:1*n]

        vrs.append(state_final[vm.indice()])
        ms.append(state_final[model.get('ina.m').indice()])
        h = state_final[model.get('ina.h').indice()]
        j = state_final[model.get('ina.j').indice()]
        avs.append(h * j)
        del(state_final)
    else:
        print(f'File not found {path}')
        vrs.append(float('nan'))
        ms.append(float('nan'))
        avs.append(float('nan'))

#
# Create figure
#
fig = plt.figure(figsize=(9, 6))  # Two-column size
fig.subplots_adjust(0.07, 0.073, 0.94, 0.95)
grid = matplotlib.gridspec.GridSpec(3, 2, wspace=0.35, hspace=0.5)

# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.98, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.004, 0.93, 'A', letter_font)
fig.text(0.004, 0.600, 'B', letter_font)
fig.text(0.470, 0.600, 'C', letter_font)
fig.text(0.470, 0.270, 'D', letter_font)

# After text
if n_beats == 1:
    atext = 'After 1 second'
else:
    atext = n_beats / 60
    atext = int(atext) if int(atext) == atext else round(atext, 1)
    atext = f'After {atext} minutes'

#
# Top: Refractory period
#
row = matplotlib.gridspec.GridSpecFromSubplotSpec(
    1, nk, subplot_spec=grid[0, :], wspace=0.2)

for i, level in enumerate(shared.ko_levels):
    print(f'Creating plot for ko={level}')
    level_str = str(float(level)).replace('.', '-')

    dname = f'{name2}-ko-{level_str}-beat{i_final}-cv.bin'
    path = os.path.join(dpath, dname)
    if not os.path.isfile(path):
        print(f'File not found: {path}')
        continue
    d = myokit.DataBlock1d.load(path).to_log().npview()

    dname = f'{name2}-ko-{level_str}-beat{i_final}-refractory-period.csv'
    path = os.path.join(dpath, dname)
    if not os.path.isfile(path):
        print(f'File not found: {path}')
        continue
    csv = myokit.DataLog.load_csv(path).npview()

    ax1 = fig.add_subplot(row[i])
    ax2 = ax1.twinx()
    ax2.set_ylim(-5, 85)

    ax1.set_xlabel('Time (s)')
    ax1.set_ylim(-95, 35)
    if i == 0:
        ax1.set_ylabel('V cell 1 (mV)')
    else:
        ax1.set_yticklabels([])
    if i == 3:
        ax2.set_ylabel('CV (cm/s)')
    else:
        ax2.set_yticklabels([' '])

    ax2.axhline(40, color='k', alpha=0.1)
    ax1.plot(d.time() * 1e-3, d['0.' + vm.qname()], color=cs[i])
    ax2.plot(csv['ts'] * 1e-3, csv['vs'], 'x-', color='#777777')

    ticks = [int(d.time()[0] * 1e-3)] * 3
    ticks[1] += cl * 1e-3 / 2
    ticks[2] += int(cl * 1e-3)
    ax1.set_xticks(ticks)
    ax1.set_xticklabels([str(x) for x in ticks])



    # Add label showing ko level
    ax2.text(0.15, 0.85, f'{level} mM', transform=ax2.transAxes)

#
# Bottom left: FRP, CV, WL
#
row = matplotlib.gridspec.GridSpecFromSubplotSpec(
    3, 1, subplot_spec=grid[1:, 0], hspace=0.1)

# Gather data
rp = []
cv = []
wl = []
sd = []
for level in ks:
    level = float(level)
    level_str = str(level).replace('.', '-')

    dname = f'{name1}-cvt-40-ko-{level_str}-beat{i_final}-wavelength.csv'
    path = os.path.join(dpath, dname)
    if os.path.isfile(path):
        print(f'Reading {path}')
        d = myokit.DataLog.load_csv(path)
        rp.append(d['rp'][0])   # ms
        cv.append(d['cv'][0])   # cm/s
        wl.append(d['wl'][0])   # mm
    else:
        print(f'File not found: {path}')
        rp.append(float('nan'))
        cv.append(float('nan'))
        wl.append(float('nan'))

    path = os.path.join(spath, f'sd-{name}-ko-{level_str}-b-{n_beats}.csv')
    if os.path.isfile(path):
        print(f'Reading {path}')
        sd.append(myokit.DataLog.load_csv(path))
    else:
        print(f'File not found: {path}')
        sd.append(float('nan'))

# Graph FRP, CV and WL
ax = fig.add_subplot(row[0, 0])
ax.set_ylabel('FRP (ms)')
ax.set_xlim(2, 8.2)
ax.set_ylim(200, 600)
ax.set_xticks(shared.ko_levels)
ax.set_xticklabels([])
ax.set_yticks([200, 300, 400, 500, 600])
ax.plot(ks, rp, '-')
for x, y, c in zip(ks, rp, shared.ko_colors):
    ax.plot(x, y, 's', color=c)
ax.text(0.98, 0.83, atext, transform=ax.transAxes, ha='right')

ax = fig.add_subplot(row[1, 0])
ax.set_ylabel('CV (cm/s)')
ax.set_xlim(2, 8.2)
ax.set_ylim(-10, 86)
ax.set_xticks(shared.ko_levels)
ax.set_xticklabels([])
ax.set_yticks([0, 20, 40, 60, 80])
ax.plot(ks, cv, '-')
for x, y, c in zip(ks, cv, shared.ko_colors):
    ax.plot(x, y, 's', color=c)

ax = fig.add_subplot(row[2, 0])
ax.set_xlabel('External potassium concentration')
ax.set_ylabel('WL (mm)')
ax.set_xlim(2, 8.2)
ax.set_ylim(160, 400)
ax.set_yticks([200, 300, 400])
ax.set_xticks(shared.ko_levels)
ax.plot(ks, wl, '-')
for x, y, c in zip(ks, wl, shared.ko_colors):
    ax.plot(x, y, 's', color=c)

#
# Bottom right: INa
#
row = matplotlib.gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=grid[1:, 1], hspace=0.5)

# INa inactivation
ax = fig.add_subplot(row[0, 0])
ax.set_ylabel('Recovered INa')
ax.set_xlabel('V (mV)')
ax.set_xlim(-130, -20)
ax.set_ylim(-0.05, 1.05)
ax.plot(vs, ss, color='#999999')
ax.plot(vs, fh(vs), '--', color='#dddddd')
#ax.plot(vs, fj(vs), ':', color='#dddddd')
for vr, k, c in zip(vrs, ks, cs):
    ax.plot([vr], [f(vr)], 'o', markersize=8, fillstyle='none', color=c,
            label=f'{k} mM')
#for av, k, c in zip(avs, ks, cs):
#    ax.axhline(av, color=c)
ax.legend(loc='lower left').get_frame().set_alpha(1)
ax.text(0.96, 0.83, atext, transform=ax.transAxes, ha='right')

# Strength-duration
ax = fig.add_subplot(row[1, 0])
ax.set_xlabel('Stimulus duration (ms)')
ax.set_ylabel('-1x Stimulus\namplitude (A/F)')
ax.set_xlim(0, 6)
ax.set_ylim(0, 42)
ax.plot(vs, ss, color='#999999')
for k, c, d in zip(ks, cs, sd):
    ax.plot(d['duration'], d['amplitude'], color=c)
ax.text(0.96, 0.83, atext, transform=ax.transAxes, ha='right')

#
# Store
#
fname = f'{fname}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')

