#!/usr/bin/env python3
#
# Figure 3: Short-term changes in membrane resistance
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
fname = 'figure-3'
debug = 'debug' in sys.argv

# Load model
model = shared.model('voigt')
name = model.name()

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 700

# Time to pre-pace after changing [K]o level
pre_time = 1000

# Prepare model
shared.prepare_model(model, protocol, fix_cleft_ko=True, pre_pace=not debug)

# Get some variables
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Perform simulations
ks = shared.ko_levels
cs = shared.ko_colors
nk = len(ks)
if debug:
    ks = [8] * nk

ds = [] # Logs with an AP in
rm = [] # Vectors time, Rm
vr = [] # Vrest, after 1 second

for k in ks:
    print(f'Simulating {k}mM')

    # Create new simulation (otherwise pre() can't be used!)
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    s.set_constant(ko, k)
    if pre_time:
        s.pre(pre_time)
    vr.append(s.state()[vm.indice()])
    d = s.run(tmax).npview()
    ds.append(d)

    dt = 0.1
    dv = 10
    rm.append(shared.rm(model, protocol, dt, dv, tmax, simulation=s))
    del(s)
    if debug:
        for i in range(1, nk):
            ds.append(ds[0])
            rm.append(rm[0])
            vr.append(vr[0])
        break


#
# Create figure
#
fig = plt.figure(figsize=(9, 4.2))  # Two-column size
fig.subplots_adjust(0.10, 0.11, 0.995, 0.935)
grid = matplotlib.gridspec.GridSpec(2, nk, wspace=0.04, hspace=0.04)

# Add model name
# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.965, shared.fancy_name(model), name_font)

letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.003, 0.880, 'A', letter_font)
fig.text(0.003, 0.470, 'B', letter_font)

def interpt(t, v, x):
    """Find time in ``t`` where ``v`` crosses point ``x``."""
    # Get last point where x is greater than x
    i = np.where(v > x)[0][-1]
    # Interpolate
    t1, t2 = t[i], t[i + 1]
    v1, v2 = v[i], v[i + 1]
    return t1 + (t2 - t1) * (x - v1) / (v2 - v1)


def interpy(t, v, tx):
    """Find value in ``v`` where ``t`` crosses point ``tx``."""
    # Get last point where t is less than tx
    i = np.where(t < tx)[0][-1]
    # Interpolate
    t1, t2 = t[i], t[i + 1]
    v1, v2 = v[i], v[i + 1]
    return v1 + (v2 - v1) * (tx - t1) / (t2 - t1)


# A: Change in AP waveform
txs = []
rxs = []
gray = '#bbbbbb'
axes = [fig.add_subplot(grid[0, i]) for i in range(nk)]
axes[0].tick_params(axis='x', direction='in')
axes[0].set_ylabel('V (mV)')
for ax, k, d, c in zip(axes, ks, ds, cs):
    t = pre_time + d.time()
    ax.plot(t, d[vm], label=str(k) + ' mM', color=c)
    ax.set_xlim(pre_time, pre_time + tmax)
    ax.set_ylim(-91, 31)
    ax.set_xticklabels('')
    ax.set_yticks([-80, -60, -40, -20, 0, 20])
    ax.legend(loc='upper right')
    t1, t2 = None, None
    try:
        t1 = interpt(t, d[vm], -20)
        ax.axvline(t1, ls='-', lw=1, color=gray, zorder=-1)
        ax.axhline(-20, ls='-', lw=1, color=gray, zorder=-1)
        t2 = interpt(t, d[vm], -40)
        ax.axvline(t2, ls='-', lw=1, color=gray, zorder=-1)
        ax.axhline(-40, ls='-', lw=1, color=gray, zorder=-1)
    except IndexError:
        pass
    txs.append((t1, t2))
    if debug:
        for i in range(1, nk):
            txs.append(txs[0])
        break

for ax in axes[1:]:
    ax.set_yticklabels('')
    ax.tick_params(axis='both', direction='in')

# B: Change in input resistance
axes = [fig.add_subplot(grid[1, i]) for i in range(nk)]
for ax in axes:
    ax.set_xlabel('Time (ms)')
    ax.set_xlim(pre_time, pre_time + tmax)
    ax.set_ylim(-250, 3200)
for ax in axes[1:]:
    ax.set_yticklabels('')
    ax.tick_params(axis='y', direction='in')
axes[0].set_ylabel('R (MOhm)')

for ax, rtr, c, tx in zip(axes, rm, cs, txs):
    rt = rtr[0] + pre_time
    rr = np.array(rtr[1])
    x = []
    for t in tx:
        if t is None:
            x.append(None)
        else:
            x.append(interpy(rt, rr, t))
    rxs.append(x)

    ax.plot(rt, rr, ':', color=c)
    rr[rr < 0] = float('nan')
    ax.plot(rt, rr, color=c)

    for t, r in zip(tx, x):
        if t is None:
            continue
        ax.axvline(t, ls='-', lw=1, color=gray, zorder=-1)
        ax.axhline(r, ls='-', lw=1, color=gray, zorder=-1)

# Align labels
fig.align_labels()

# Annotate
ar = {'arrowstyle': '->'}

# Annotate graph for 2.5
level = 2.5
if level in ks:
    i = ks.index(level)
    #t, r = 1010, rm[i][1][0]
    #axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t + 5, r + 1000), arrowprops=ar)
    t, r = txs[i][0], rxs[i][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 30, r + 700), arrowprops=ar)
    t, r = txs[i][1], rxs[i][1]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t + 10, r + 1000), arrowprops=ar)

# Annotate graph for 3.2
level = 3.2
if level in ks:
    i = ks.index(level)
    #t, r = 1010, rm[i][1][0]
    #axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t + 5, r + 1000), arrowprops=ar)
    t, r = txs[i][0], rxs[i][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 60, r + 700), arrowprops=ar)
    t, r = txs[i][1], rxs[i][1]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 10, r + 1000), arrowprops=ar)

# Annotate graph for 4
level = 4
if level in ks:
    i = ks.index(level)
    t, r = 1010, rm[i][1][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t + 5, r + 1000), arrowprops=ar)
    t, r = txs[i][0], rxs[i][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 30, r + 700), arrowprops=ar)
    t, r = txs[i][1], rxs[i][1]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 40, r + 900), arrowprops=ar)

# Annotate graph for 5.4
level = 5.4
if level in ks:
    i = ks.index(level)
    t, r = 1010, rm[i][1][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t + 5, r + 1100), arrowprops=ar)
    t, r = txs[i][0], rxs[i][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 30, r + 700), arrowprops=ar)
    t, r = txs[i][1], rxs[i][1]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t + 10, r + 1200), arrowprops=ar)

# Annotate graph for 8
level = 8
if level in ks:
    i = ks.index(level)
    t, r = 1010, rm[i][1][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t, r + 1200), arrowprops=ar)
    t, r = txs[i][0], rxs[i][0]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r), (t - 40, r + 700), arrowprops=ar)
    t, r = txs[i][1], rxs[i][1]
    axes[i].annotate(str(int(np.round(r, 0))), (t, r + 50), (t + 30, r + 1000), arrowprops=ar)


# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
