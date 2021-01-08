#!/usr/bin/env python3
#
# Model figure m5: Strand simulation results
#
import os
import matplotlib.pyplot as plt
from  matplotlib.patches import Rectangle
import myokit
import shared

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-m6c-strand'
#debug = 'debug' in sys.argv

# Get path for data
dpath = 'strand-data'

# Load model, just for name
model = shared.model()
name = model.name()

# Gather data
ko = []
rp = []
cv = []
wl = []
for level in shared.ko_levels:
    level = float(level)
    level_str = str(level).replace('.', '-')

    dname = f'{name}-cvt-40-ko-{level_str}-beat52-wavelength.csv'
    path = os.path.join(dpath, dname)
    if os.path.isfile(path):
        print('Reading: ' + path)
        d = myokit.DataLog.load_csv(path)
        ko.append(level)
        rp.append(d['rp'][0])   # ms
        cv.append(d['cv'][0])   # cm/s
        wl.append(d['wl'][0])   # mm
    else:
        print('File not found: ' + path)

#
# Figure without AF
#
xlim = min(shared.ko_levels) - 0.2, max(shared.ko_levels) + 0.2

fig = plt.figure(figsize=(4.3, 4))  # One-column size
fig.subplots_adjust(0.19, 0.12, 0.91, 0.98, 0, 0.1)

# Add model name
fig.patches.append(Rectangle((0.93, 0), 0.07, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.95, 0.5, name, name_font)

# Graph AP and CV
ax1 = fig.add_subplot(3, 1, 1)
ax1.set_ylabel('FRP (ms)')
ax1.set_xlim(*xlim)
ax1.set_xticks(shared.ko_levels)
ax1.set_xticklabels([])
ax1.plot(ko, rp, '-')
for x, y, c in zip(ko, rp, shared.ko_colors):
    ax1.plot(x, y, 's', color=c)

ax2 = fig.add_subplot(3, 1, 2)
ax2.set_ylabel('CV (cm/s)')
ax2.set_xlim(*xlim)
ax2.set_xticks(shared.ko_levels)
ax2.set_xticklabels([])
ax2.plot(ko, cv, '-')
for x, y, c in zip(ko, cv, shared.ko_colors):
    ax2.plot(x, y, 's', color=c)

ax3 = fig.add_subplot(3, 1, 3)
ax3.set_xlabel('External potassium')
ax3.set_ylabel('WL (mm)')
ax3.set_xlim(*xlim)
ax3.set_xticks(shared.ko_levels)
ax3.plot(ko, wl, '-')
for x, y, c in zip(ko, wl, shared.ko_colors):
    ax3.plot(x, y, 's', color=c)

path = os.path.join(fpath, fname + '-' + name)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')



#
# Figure with AF
#
if False:
    # Gather data
    ko_af = []
    rp_af = []
    cv_af = []
    wl_af = []

    for level in shared.ko_levels:
        level = float(level)
        level_str = str(level).replace('.', '-')

        dname = f'{name}-af-cvt-40-ko-{level_str}-wavelength.csv'
        path = os.path.join(dpath, dname)
        if os.path.isfile(path):
            print('Reading: ' + path)
            d = myokit.DataLog.load_csv(path)
            ko_af.append(level)
            rp_af.append(d['rp'][0])   # ms
            cv_af.append(d['cv'][0])   # cm/s
            wl_af.append(d['wl'][0])   # mm
        else:
            print('File not found: ' + path)


    fig = plt.figure(figsize=(4.3, 4))  # One-column size
    fig.subplots_adjust(0.19, 0.12, 0.91, 0.92, 0, 0.1)

    # Add model name
    fig.patches.append(Rectangle((0.93, 0), 0.07, 1, color='#cccccc', fill=True,
        transform=fig.transFigure, zorder=-1))
    name_font = {
        'fontsize': 12,
        'rotation': 'vertical',
        'verticalalignment': 'center',
    }
    fig.text(0.95, 0.5, name, name_font)

    # Graph AP and CV
    ax1 = fig.add_subplot(3, 1, 1)
    ax1.set_ylabel('FRP (ms)')
    ax1.set_xlim(*xlim)
    ax1.set_xticks(shared.ko_levels)
    ax1.set_xticklabels([])
    ax1.plot(ko, rp, 's-', label='baseline')
    ax1.plot(ko_af, rp_af, 's-', label='AF')
    ax1.legend(ncol=2, loc=(0.2, 1.05))

    ax2 = fig.add_subplot(3, 1, 2)
    ax2.set_ylabel('CV (cm/s)')
    ax2.set_xlim(*xlim)
    ax2.set_xticks(shared.ko_levels)
    ax2.set_xticklabels([])
    ax2.plot(ko, cv, 's-')
    ax2.plot(ko_af, cv_af, 's-')

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.set_xlabel('External potassium')
    ax3.set_ylabel('WL (mm)')
    ax3.set_xlim(*xlim)
    ax3.set_xticks(shared.ko_levels)
    ax3.plot(ko, wl, 's-')
    ax3.plot(ko_af, wl_af, 's-')

    path = os.path.join(fpath, fname + '-' + name + '-af')
    plt.savefig(path + '.png')
    plt.savefig(path + '.pdf')

