#!/usr/bin/env python3
#
# Combine strand simulation results for multiple concentrations
#
import os
import matplotlib
import matplotlib.pyplot as plt
from  matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import myokit
import shared

matplotlib.rcParams['axes.spines.right'] = True

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-m6b-frp'
#debug = 'debug' in sys.argv

# Settings
af = False
cvt = 40

# Load model
model = shared.model()
print(f'Loaded {model.name()}')

# Get data file names
dpath = 'strand-data'
name1 = model.name() + ('-af' if af else '')
name2 = f'{name1}-cvt-{cvt}'

# Set af mode
if af:
    print('Enabling AF mode')
    shared.enable_af(model)

# Get some variables
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Ko levels and colours
ks = shared.ko_levels
cs = shared.ko_colors
n = len(ks)

# Create figure
fig = plt.figure(figsize=(9, 2.3))  # Two-column size
fig.subplots_adjust(0.07, 0.20, 0.91, 0.98)
grid = GridSpec(1, n, wspace=0.2)

# Add model name
fig.patches.append(Rectangle((0.965, 0), 0.035, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.975, 0.5, name1, name_font)

# Add panels
for i, level in enumerate(shared.ko_levels):
    print(f'Creating plot for ko={level}')
    level_str = str(float(level)).replace('.', '-')

    dname = f'{name2}-ko-{level_str}-beat52-cv.bin'
    path = os.path.join(dpath, dname)
    d = myokit.DataBlock1d.load(path).to_log()

    dname = f'{name2}-ko-{level_str}-beat52-refractory-period.csv'
    path = os.path.join(dpath, dname)
    csv = myokit.DataLog.load_csv(path)

    ax1 = plt.subplot(grid[i])
    ax2 = ax1.twinx()
    ax2.set_ylim(-5, 85)

    ax1.set_xlabel('Time (ms)')
    ax1.set_ylim(-95, 45)
    if i == 0:
        ax1.set_ylabel('Vm in cell 1 (mV)')
    else:
        ax1.set_yticklabels([])
    if i == n - 1:
        ax2.set_ylabel('CV (cm/s)')
    else:
        ax2.set_yticklabels([])

    ax2.axhline(40, color='k', alpha=0.1)
    ax1.plot(d.time(), d['0.' + vm.qname()], color=cs[i])
    ax2.plot(csv['ts'], csv['vs'], 'x-', color='#777777')

    # Add label showing ko level
    #ypos = 0.55 if csv['vs'][-1] > 60 else 0.9
    ypos = 0.40
    ax2.text(0.66, ypos, f'{level} mM', transform=ax2.transAxes)

# Store
fname = f'{fname}-{name1}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')

