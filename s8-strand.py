#!/usr/bin/env python3
#
# Strand simulation figures
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import os
import shared
import sys

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-s8-strand'

# Load model
model = shared.model()
shared.splash(fname, model.name())
# Note: Not preparing model here, that will be done in strand sim call

# Perform test beat if asked
if '-g' in sys.argv:
    shared.strand_test_beat(model, conductance=True)
    sys.exit(0)
elif '-f' in sys.argv:
    shared.strand_test_beat(model, conductance=False)
    sys.exit(0)
elif '-s' in sys.argv:
    shared.strand_test_stability(model)
    sys.exit(0)

# Show crosses instead of squares at Ko levels without AP
if 'voigt' in model.name():
    crosses = [2.5]
elif 'grandi' in model.name():
    crosses = [2.5, 3.2]
elif 'nygren' in model.name():
    crosses = [2.5]
elif 'koivumaki' in model.name():
    crosses = [2.5, 3.2]
else:
    crosses = []

# Create figure
fig = plt.figure(figsize=(9, 4.5))  # Two-column size
fig.subplots_adjust(0.08, 0.11, 0.95, 0.98)
grid = matplotlib.gridspec.GridSpec(1, 2, wspace=0.35, hspace=0.3)

# Add model name
fig.patches.append(matplotlib.patches.Rectangle(
    (0.965, 0), 0.035, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.975, 0.5, model.name(), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}

# Top: Refractory period
axa = fig.add_subplot(grid[0, 0])
shared.strand_frp_graph(axa, model)

# Bottom left: FRP, CV, WL
hs = 0.3
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    5, 1, subplot_spec=grid[0, 1], hspace=hs)
ax1 = fig.add_subplot(sg[0, 0])
ax2 = fig.add_subplot(sg[1, 0])
ax3 = fig.add_subplot(sg[2, 0])
ax4 = fig.add_subplot(sg[3, 0])
ax5 = fig.add_subplot(sg[4, 0])
shared.strand_stats_graphs(ax1, ax2, ax3, ax4, ax5, model, False, crosses)
fig.align_ylabels([ax1, ax2, ax3, ax4, ax5])

# Store
fname = f'{fname}-{model.name()}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')
print('Done')
