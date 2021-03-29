#!/usr/bin/env python3
#
# Figure 5: Combine strand simulation results for multiple concentrations
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import os
import shared

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-5'
shared.splash(fname)

# Load model
model = shared.model('voigt')
print(f'Loaded {model.name()}')
# Note: Not preparing model here, that will be done in strand sim call

# Create figure
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

# Top: Refractory period
axa = fig.add_subplot(grid[0, 0])
shared.strand_frp_graph(axa, model, False)
axa.text(0.5, 1.05, 'After 1 second', transform=axa.transAxes, ha='center')

axb = fig.add_subplot(grid[0, 1])
shared.strand_frp_graph(axb, model, True)
axb.text(0.5, 1.05, 'After 15 minutes', transform=axb.transAxes, ha='center')


# Bottom left: FRP, CV, WL
hs = 0.3
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    5, 1, subplot_spec=grid[1:, 0], hspace=hs)
ax1 = fig.add_subplot(sg[0, 0])
ax2 = fig.add_subplot(sg[1, 0])
ax3 = fig.add_subplot(sg[2, 0])
ax4 = fig.add_subplot(sg[3, 0])
ax5 = fig.add_subplot(sg[4, 0])
shared.strand_stats_graphs(
    ax1, ax2, ax3, ax4, ax5, model, False)
fig.align_ylabels([axa, ax1, ax2, ax3, ax4, ax5])

sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    5, 1, subplot_spec=grid[1:, 1], hspace=hs)
ax1 = fig.add_subplot(sg[0, 0])
ax2 = fig.add_subplot(sg[1, 0])
ax3 = fig.add_subplot(sg[2, 0])
ax4 = fig.add_subplot(sg[3, 0])
ax5 = fig.add_subplot(sg[4, 0])
shared.strand_stats_graphs(
    ax1, ax2, ax3, ax4, ax5, model, True)
fig.align_ylabels([axb, ax1, ax2, ax3, ax4, ax5])

# Store
fname = f'{fname}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')
plt.savefig(path + '.jpg', dpi=600)
print('Done')
