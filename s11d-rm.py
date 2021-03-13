#!/usr/bin/env python3
#
# Short-term changes in membrane resistance
#
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import matplotlib.patches
import matplotlib.transforms
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s11d-rm'

# Load model
model = shared.model()
name = model.name()
shared.splash(fname, name)

# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Maximum time to show in plots
tmax = 700

# Prepare model
shared.prepare_model(model, protocol)

# Get some variables
tt = model.time()
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Get concentrations
ks = shared.ko_levels
cs = shared.ko_colors
nk = len(ks)

# Gather data
dt = 0.1
dv = 10
ds = []  # Logs with an AP in
rs = []  # Logs with 'time' and 'resistance'

for k in ks:
    dname = os.path.join(shared.sdata, f'rm-{name}')
    dname += f'-dt-{shared.fstr(dt)}'
    dname += f'-dv-{shared.fstr(dv)}'
    dname += f'-ko-{shared.fstr(k)}'
    path1 = f'{dname}-ap.zip'
    path2 = f'{dname}-rm.zip'
    if os.path.isfile(path1) and os.path.isfile(path2):
        print(f'Loading cached data for {k}')
        ds.append(myokit.DataLog.load(path1).npview())
        rs.append(myokit.DataLog.load(path2).npview())
    else:
        print(f'Simulating {k}mM')

        # Create new simulation (otherwise pre() can't be used!)
        s = myokit.Simulation(model, protocol)
        s.set_tolerance(1e-8, 1e-8)
        s.set_constant(ko, k)
        s.pre(cl)
        d = s.run(tmax, log=[tt, vm]).npview()
        d[tt] += cl
        ds.append(d)
        d.save(path1)

        r = shared.rm(model, protocol, dt, dv, tmax, simulation=s)
        d = myokit.DataLog(time='time')
        d['time'] = r[0] + cl
        d['resistance'] = r[1]
        rs.append(d)
        d.save(path2)
        del(s, r, d)

vr = [d[vm][0] for d in ds]


#
# Create figure
#
fig = plt.figure(figsize=(9, 4.2))  # Two-column size
fig.subplots_adjust(0.08, 0.11, 0.95, 0.98)
grid = matplotlib.gridspec.GridSpec(2, nk, wspace=0.08, hspace=0.05)

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


def interpt(t, v, x):
    """Find time in ``t`` where ``v`` crosses point ``x``."""
    # Get last point where x is greater than x
    i = np.where(v > x)[0][-1]
    # Edge case
    if i + 1 == len(t):
        return t[-1]
    # Interpolate
    t1, t2 = t[i], t[i + 1]
    v1, v2 = v[i], v[i + 1]
    return t1 + (t2 - t1) * (x - v1) / (v2 - v1)


def interpy(t, v, tx):
    """Find value in ``v`` where ``t`` crosses point ``tx``."""
    # Get last point where t is less than tx
    i = np.where(t < tx)[0][-1]
    # Edge case
    if i + 1 == len(t):
        return v[-1]
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
    ax.plot(d[tt], d[vm], label=str(k) + ' mM', color=c)
    ax.set_xlim(cl, cl + tmax)
    ax.set_xticklabels('')
    ax.set_ylim(-105, 31)
    #ax.set_ylim(-81, 21)
    #ax.set_yticks([-80, -60, -40, -20, 0, 20])
    ax.legend(loc='upper right')
    t1, t2 = None, None
    try:
        t1 = interpt(d[tt], d[vm], -20)
        ax.axvline(t1, ls='-', lw=1, color=gray, zorder=-1)
        ax.axhline(-20, ls='-', lw=1, color=gray, zorder=-1)
        t2 = interpt(d[tt], d[vm], -40)
        ax.axvline(t2, ls='-', lw=1, color=gray, zorder=-1)
        ax.axhline(-40, ls='-', lw=1, color=gray, zorder=-1)
    except IndexError:
        pass
    txs.append((t1, t2))

for ax in axes[1:]:
    ax.set_yticklabels('')
    #ax.tick_params(axis='both', direction='in')

# B: Change in input resistance
axes = [fig.add_subplot(grid[1, i]) for i in range(nk)]
for ax in axes:
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_xlabel('Time (ms)')
    ax.set_xlim(cl, cl + tmax)
    if 'voigt' in model.name():
        ax.set_ylim(-250, 3200)
    else:
        ax.set_ylim(-250, 3200)

for ax in axes[1:]:
    ax.set_yticklabels('')
    #ax.tick_params(axis='y', which='both', direction='in')
axes[0].set_ylabel('R (M$\Omega$)')

# For each concentration...
for ax, d, c, tx in zip(axes, rs, cs, txs):
    rt = d.time()
    rr = np.array(d['resistance'])
    x = []
    for t in tx:
        if t is None:
            x.append(None)
        else:
            x.append(interpy(rt, rr, t))
    rxs.append(x)

    # Plot straight line at any zero crossing
    ic = np.nonzero(rr[1:] * rr[:-1] < 0)
    for i in ic:
        ax.plot((rt[i + 1], rt[i + 1]), (rr[i], rr[i + 1]), ':', color=c)

    # Shade areas with negative resistance
    trans = matplotlib.transforms.blended_transform_factory(
        ax.transData, ax.transAxes)
    ax.fill_between(rt, 0, 1, where=rr <= 0,
                    facecolor='black', alpha=0.1, transform=trans)

    # Plot line with zero parts filtered out
    #rr[rr < 0] = float('nan')
    ax.plot(rt, rr, color=c)

    for t, r in zip(tx, x):
        if t is None:
            continue
        ax.axvline(t, ls='-', lw=1, color=gray, zorder=-1)
        ax.axhline(r, ls='-', lw=1, color=gray, zorder=-1)


# Align labels
fig.align_labels()

# Annotate Voigt-Heijman model graph
if 'voigt' in name:
    ar = {'arrowstyle': '->'}

    # Annotate graph for 2.5
    level = 2.5
    if level in ks:
        i = ks.index(level)
        t, r = txs[i][0], rxs[i][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 30, r + 700), arrowprops=ar)
        t, r = txs[i][1], rxs[i][1]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t + 10, r + 1000),
            arrowprops=ar)

    # Annotate graph for 3.2
    level = 3.2
    if level in ks:
        i = ks.index(level)
        t, r = txs[i][0], rxs[i][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 60, r + 700), arrowprops=ar)
        t, r = txs[i][1], rxs[i][1]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 10, r + 1000),
            arrowprops=ar)

    # Annotate graph for 4
    level = 4
    if level in ks:
        i = ks.index(level)
        t, r = 1010, rs[i]['resistance'][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t + 15, r + 1000),
            arrowprops=ar)
        t, r = txs[i][0], rxs[i][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 30, r + 700), arrowprops=ar)
        t, r = txs[i][1], rxs[i][1]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 40, r + 900), arrowprops=ar)

    # Annotate graph for 5.4
    level = 5.4
    if level in ks:
        i = ks.index(level)
        t, r = 1010, rs[i]['resistance'][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t + 5, r + 1100), arrowprops=ar)
        t, r = txs[i][0], rxs[i][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 30, r + 700), arrowprops=ar)
        t, r = txs[i][1], rxs[i][1]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 10, r + 1200),
            arrowprops=ar)

    # Annotate graph for 8
    level = 8
    if level in ks:
        i = ks.index(level)
        t, r = 1010, rs[i]['resistance'][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t, r + 1200), arrowprops=ar)
        t, r = txs[i][0], rxs[i][0]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r), (t - 40, r + 700), arrowprops=ar)
        t, r = txs[i][1], rxs[i][1]
        axes[i].annotate(
            str(int(np.round(r, 0))), (t, r + 50), (t + 30, r + 1000),
            arrowprops=ar)


# Show / store
path = os.path.join(fpath, f'{fname}-{name}')
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
