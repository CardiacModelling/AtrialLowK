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
fname = 'figure-7'
debug = 'debug' in sys.argv

# Settings
cvt = 40
i1 = 52
i2 = 951

# Load model
model = shared.model('voigt')
name = model.name()
print(f'Loaded {name}')

# Get data file names
dpath = 'strand-data'
spath = 'sd-data'
name2 = f'{name}-cvt-{cvt}'

# Create protocol
cl = 1000
offset = 50

# Prepare model
shared.prepare_model(model, fix_cleft_ko=True, pre_pace=False)

# Get some variables
vm = model.labelx('membrane_potential')

# Ko levels and colours
ks = shared.ko_levels
cs = shared.ko_colors
nk = len(ks)


#
# Create figure
#
fig = plt.figure(figsize=(9, 8))  # Two-column size
fig.subplots_adjust(0.10, 0.060, 0.98, 0.92)
grid = matplotlib.gridspec.GridSpec(3, 2, wspace=0.35, hspace=0.3)

# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.98, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.003, 0.940, 'A', letter_font)
fig.text(0.003, 0.620, 'B', letter_font)
fig.text(0.003, 0.500, 'C', letter_font)
fig.text(0.003, 0.380, 'D', letter_font)
fig.text(0.003, 0.260, 'E', letter_font)
fig.text(0.003, 0.140, 'F', letter_font)

#
# Top: Refractory period
#
def refrac(ax, i_final, atext):
    # Plot refractory period exploration graphs

    # Gather data
    ap, rp, cv, wl, vr = [], [], [], [], []

    # Update axes
    ax.minorticks_on()
    ax.set_xlabel('Time (s)')
    ax.set_ylim(-85, 35)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('V cell 1 (mV)')

    for i, level in enumerate(shared.ko_levels):
        print(f'Creating plot for ko={level}')
        level_str = str(float(level)).replace('.', '-')

        # Load AP data
        dname = f'{name2}-ko-{level_str}-beat{i_final}-cv.bin'
        path = os.path.join(dpath, dname)
        if not os.path.isfile(path):
            print(f'File not found: {path}')
            continue
        d = myokit.DataBlock1d.load(path).to_log().npview()

        # Gather data
        dname = f'{name}-cvt-40-ko-{level_str}-beat{i_final}-wavelength.csv'
        path = os.path.join(dpath, dname)
        if os.path.isfile(path):
            print(f'Reading {path}')
            csv = myokit.DataLog.load_csv(path)
            ap.append(csv['apd'][0])  # ms
            rp.append(csv['rp'][0])   # ms
            cv.append(csv['cv'][0])   # cm/s
            wl.append(csv['wl'][0])   # mm
        else:
            print(f'File not found: {path}')
            ap.append(float('nan'))
            rp.append(float('nan'))
            cv.append(float('nan'))
            wl.append(float('nan'))

        # Plot AP
        ax.plot(d.time() * 1e-3, d[vm, 0], color=cs[i], label=f'{level} mM')

        # Plot FRP
        if np.isfinite(rp[-1]):
            frp = (d.time()[0] + offset + rp[-1]) * 1e-3
            ax.axvline(frp, color=cs[i], alpha=0.5)

        # Store vrest
        vr.append(d[vm, 0][0])

    # Finalise plot
    ticks = [int(d.time()[0] * 1e-3)] * 3
    ticks[1] += cl * 1e-3 / 2
    ticks[2] += int(cl * 1e-3)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) for x in ticks])
    ax.legend(loc='upper right')

    ax.text(0.5, 1.05, atext, transform=ax.transAxes, ha='center')

    return ap, rp, cv, wl, vr

axa = fig.add_subplot(grid[0, 0])
data1 = refrac(axa, i1, 'After 1 second')

axb = fig.add_subplot(grid[0, 1])
data2 = refrac(axb, i2, 'After 15 minutes')

#
# Bottom left: FRP, CV, WL
#
def summary(ax1, ax2, ax3, ax4, ax5, data):

    # Unpack
    ap, rp, cv, wl, vr = data

    # Shared limits
    xlim = 2.3, 5.6

    # Vr
    ax = ax1
    ax.set_ylabel('$V_r$ (mV)')
    ax.set_xlim(*xlim)
    ax.set_ylim(-78, -49)
    ax.set_xticks(shared.ko_levels)
    ax.set_xticklabels([])
    ax.plot(ks, vr, '-')
    for x, y, c in zip(ks, vr, shared.ko_colors):
        ax.plot(x, y, 's', color=c)

    # CV
    ax = ax2
    ax.set_ylabel('CV (cm/s)')
    ax.set_xlim(*xlim)
    ax.set_ylim(-10, 86)
    ax.set_xticks(shared.ko_levels)
    ax.set_xticklabels([])
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.plot(ks, cv, '-')
    for x, y, c in zip(ks, cv, shared.ko_colors):
        ax.plot(x, y, 's', color=c)

    # APD
    ax = ax3
    ax.set_ylabel('APD (ms)')
    ax.set_xlim(*xlim)
    ax.set_ylim(200, 600)
    ax.set_xticks(shared.ko_levels)
    ax.set_xticklabels([])
    ax.set_yticks([200, 300, 400, 500, 600])
    ax.plot(ks, rp, '-')
    for x, y, c in zip(ks, rp, shared.ko_colors):
        ax.plot(x, y, 's', color=c)

    # FRP
    ax = ax4
    ax.set_ylabel('FRP (ms)')
    ax.set_xlim(*xlim)
    ax.set_ylim(200, 600)
    ax.set_xticks(shared.ko_levels)
    ax.set_xticklabels([])
    ax.set_yticks([200, 300, 400, 500, 600])
    ax.plot(ks, rp, '-')
    for x, y, c in zip(ks, rp, shared.ko_colors):
        ax.plot(x, y, 's', color=c)

    # WL
    ax = ax5
    ax.set_xlabel('External potassium concentration $[K^+]_o$')
    ax.set_ylabel('WL (mm)')
    ax.set_xlim(*xlim)
    ax.set_ylim(160, 400)
    ax.set_yticks([200, 300, 400])
    ax.set_xticks(shared.ko_levels)
    ax.plot(ks, wl, '-')
    for x, y, c in zip(ks, wl, shared.ko_colors):
        ax.plot(x, y, 's', color=c)

hs = 0.3
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    5, 1, subplot_spec=grid[1:, 0], hspace=hs)
ax1 = fig.add_subplot(sg[0, 0])
ax2 = fig.add_subplot(sg[1, 0])
ax3 = fig.add_subplot(sg[2, 0])
ax4 = fig.add_subplot(sg[3, 0])
ax5 = fig.add_subplot(sg[4, 0])
summary(ax1, ax2, ax3, ax4, ax5, data1)
fig.align_ylabels([axa, ax1, ax2, ax3, ax4, ax5])

sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    5, 1, subplot_spec=grid[1:, 1], hspace=hs)
ax1 = fig.add_subplot(sg[0, 0])
ax2 = fig.add_subplot(sg[1, 0])
ax3 = fig.add_subplot(sg[2, 0])
ax4 = fig.add_subplot(sg[3, 0])
ax5 = fig.add_subplot(sg[4, 0])
summary(ax1, ax2, ax3, ax4, ax5, data2)
fig.align_ylabels([axb, ax1, ax2, ax3, ax4, ax5])

#
# Store
#
fname = f'{fname}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')

