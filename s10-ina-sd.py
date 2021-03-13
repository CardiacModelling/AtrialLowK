#!/usr/bin/env python3
#
# INa availability and SD curves
#
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.gridspec
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s10-ina-sd'
debug = 'debug' in sys.argv

# Load model
model = shared.model()
name = model.name()
shared.splash(fname, name)
shared.prepare_model(model)

# Load data
ks = shared.ko_levels
cs = shared.ko_colors
vs, ss, v1s, a1s, v2s, a2s = shared.ina_availability(model)

# Start testing
sd1, sd2 = [], []
for k in ks:
    print(f'Testing {k} mM')
    sd1.append(shared.sd_log_fast(model, k))
    sd2.append(shared.sd_log_slow(model, k))
print('Done')


#
# Create figure
#
fig = plt.figure(figsize=(9, 4.5))  # Two-column size
fig.subplots_adjust(0.07, 0.1, 0.95, 0.98)
grid = matplotlib.gridspec.GridSpec(2, 3, wspace=0.6, hspace=0.25)

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

# Top-left: INa 1 beat
ax00 = fig.add_subplot(grid[0, 0])
#ax00.set_xlabel('Voltage (mV)')
ax00.set_ylabel('Recovered INa')
ax00.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax00.set_xlim(-130, -20)
ax00.set_ylim(-0.05, 1.25)
ax00.plot(vs, ss, color='#999999')
for vr, a, k, c in zip(v1s, a1s, ks, cs):
    ax00.plot([vr], [a], 'o', markersize=8, fillstyle='none', color=c,
              label=f'{k} mM')
ax00.legend(loc=(0.65, 0.35), handlelength=1, handletextpad=0.5)
ax00.text(0.043, 0.88, 'After 1 second', transform=ax00.transAxes)

# Bottom-left: INa 900 beats
ax10 = fig.add_subplot(grid[1, 0])
ax10.set_xlabel('Voltage (mV)')
ax10.set_ylabel('Recovered INa')
ax10.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax10.set_xlim(-130, -20)
ax10.set_ylim(-0.05, 1.25)
ax10.plot(vs, ss, color='#999999')
for vr, a, k, c in zip(v2s, a2s, ks, cs):
    ax10.plot([vr], [a], 'o', markersize=8, fillstyle='none', color=c,
              label=f'{k}')
ax10.text(0.043, 0.88, 'After 15 minutes', transform=ax10.transAxes)
fig.align_labels((ax00, ax10))

# Top-right: SD curve 1 beat
ax = fig.add_subplot(grid[0, 1:])
ax.set_ylabel('-1x stimulus\namplitude (A/F)')
ax.set_xlim(0, 5)
ax.set_ylim(0, 50)
for k, c, d in zip(ks, cs, sd1):
    if d is not None:
        ls = '--' if k in [5.4] else '-'
        ax.plot(d['duration'], d['amplitude'], color=c, ls=ls, label=f'{k} mM')
#ax.text(0.96, 0.83, atext, transform=ax.transAxes, ha='right')
ax.legend(loc='upper right')
ax.text(0.10, 1, 'After 1 second')

# Bottom-right: SD curve 900 beats
ax = fig.add_subplot(grid[1, 1:])
ax.set_xlabel('Stimulus duration (ms)')
ax.set_ylabel('-1x stimulus\namplitude (A/F)')
ax.set_xlim(0, 5)
ax.set_ylim(0, 50)
for k, c, d in zip(ks, cs, sd2):
    if d is not None:
        ls = '--' if k in [5.4] else '-'
        ax.plot(d['duration'], d['amplitude'], color=c, ls=ls)
ax.text(0.10, 1, 'After 15 minutes')

# Show / store
path = os.path.join(fpath, f'{fname}-{name}')
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')

