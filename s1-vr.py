#!/usr/bin/env python3
#
# Show literature data Vr versus Ko
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import os
import shared
import sys

# Get path for figures
sim = 'sim' in sys.argv[1:]
fpath = shared.figure_dir()
if sim:
    fname = 'figure-s10-vr'
else:
    fname = 'figure-s1-vr'
shared.splash(fname)

# Consistent colors for models
cmap = matplotlib.cm.get_cmap('tab10')

# Load and adjust literature data
dpath = 'literature-data'
d1 = os.path.join(dpath, 'gelband-1972-fig-6.csv')
d1 = myokit.DataLog.load_csv(d1).npview()
d2 = os.path.join(dpath, 'ten-eick-1979-fig-4.csv')
d2 = myokit.DataLog.load_csv(d2).npview()

# Get model data
levels = [
    1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875,
    2, 2.125, 2.25, 2.375, 2.5, 2.75,
    3, 3.25, 3.5, 3.75,
    4, 4.25, 4.5, 4.75,
    5, 5.4,
    6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5,
    10, 11, 12, 12.5, 13, 14, 15,
    20, 21, 22, 23, 24, 25,
    30, 32.5, 35, 37.5, 40, 45, 50,
]

models = shared.models()
v1s = {}
v2s = {}
protocol = shared.default_protocol()
for model in models.values():
    print(f'Gathering data for {model.name()}')
    shared.prepare_model(model)
    vm = model.label('membrane_potential').indice()
    v1s[shared.fancy_name(model)] = [
        shared.state_fast(model, protocol, k, lw=True)[vm] for k in levels]
    v2s[shared.fancy_name(model)] = [
        shared.state_slow(model, protocol, k, lw=True)[vm] for k in levels]

if True:
    print(f'Gathering data for LR1-IK1 Courtemanche')
    model = myokit.load_model('models/courtemanche-lr1-ik1.mmt')
    shared.prepare_model(model)
    vm = model.label('membrane_potential').indice()
    v1s['Courtemanche with LR1 IK1'] = [
        shared.state_fast(model, protocol, k, lw=True)[vm] for k in levels]
    v2s['Courtemanche with LR1 IK1'] = [
        shared.state_slow(model, protocol, k, lw=True)[vm] for k in levels]



#
# Create figure
#
fig = plt.figure(figsize=(9, 10))  # Two-column size
#fig.subplots_adjust(0.15, 0.15, 0.98, 0.98)
fig.subplots_adjust(0.085, 0.05, 0.99, 0.995)
grid = matplotlib.gridspec.GridSpec(2, 1, hspace=0.15)

# Minor ticks
#matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True

ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel('[K$^+$]$_o$ (mM)')
ax.set_ylabel('$V_r$ (mV)')
ax.set_xscale('log')
ax.errorbar(
    d1['Ko'], d1['Vr'], d1['err'], fmt='s-', label='Gelband et al. 1972 (n=6)',
    capsize=5, zorder=98, color=cmap(8))
ax.plot(
    d2['Ko'], d2['Vr'], '^-', label='Ten Eick et al. 1979 (n=2)',
    zorder=99, color=cmap(9))
if sim:
    for name, vr in v1s.items():
        ax.plot(levels, vr, '.--', alpha=0.5, label=name)
ax.set_xticks(d1['Ko'])
ax.set_xticklabels([int(x) for x in d1['Ko']])
ax.legend(loc='lower right', ncol=2)
if sim:
    ax.text(0.02, 0.95, 'Simulations 1 second after K+ change',
            transform=ax.transAxes)

ax = fig.add_subplot(grid[1, 0])
ax.set_xlabel('[K$^+$]$_o$ (mM)')
ax.set_ylabel('$V_r$ (mV)')
if sim:
    ax.set_xscale('log')
ax.errorbar(
    d1['Ko'], d1['Vr'], d1['err'], fmt='s-', label='Gelband et al. 1972 (n=6)',
    capsize=5, zorder=98, color=cmap(8))
ax.plot(
    d2['Ko'], d2['Vr'], '^-', label='Ten Eick et al. 1979 (n=2)',
    zorder=99, color=cmap(9))
if sim:
    for name, vr in v2s.items():
        ax.plot(levels, vr, '.--', alpha=0.5, label=name)
ax.set_xticks(d1['Ko'])
ax.set_xticklabels([int(x) for x in d1['Ko']])
ax.legend(loc='lower right', ncol=2)
if sim:
    ax.text(0.02, 0.95, 'Simulations 15 minutes after K+ change',
            transform=ax.transAxes)


#
# Store
#
fname = f'{fname}'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.png')
plt.savefig(f'{path}.pdf')

