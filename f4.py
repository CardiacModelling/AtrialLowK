#!/usr/bin/env python3
#
# Figure 4: Detailed changes in voigt after 8 minutes
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import shared

# Minor y-ticks all round
matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-4'
debug = 'debug' in sys.argv
shared.splash(fname)

# Load model
model = shared.model('voigt')
name = model.name()

# Create protocol
cl = 1000   # MUST be 1000 for units nC to work out
protocol = shared.default_protocol()

# Time to pre-pace after changing [K]o level
n_beats = 481
text_1 = '1 sec'
text_2 = '8 min'
pre_time_1 = cl
pre_time_2 = n_beats * cl

# Prepare model
shared.prepare_model(model, protocol, pA=True)

# Get some variables
time = model.time()
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')
nai = model.labelx('Na_i')
ina = model.labelx('I_Na')
inal = model.labelx('I_NaL')
inab = model.labelx('I_NaB')
inak = model.labelx('I_NaK')
inaca = model.labelx('I_NaCa')

# Get steady-state inactivation curve
ina_j_inf = model.labelx('ina_j_inf')
ina_h_inf = model.labelx('ina_h_inf')
fj = ina_j_inf.pyfunc()
fh = ina_h_inf.pyfunc()
f = lambda v: fj(v) * fh(v)
ftm = model.get('ina.m.tau').pyfunc()
vs = np.linspace(-160, 40, 200)
ss = f(vs)

# Perform simulations
s = myokit.Simulation(model, protocol)
ks = shared.ko_levels
cs = shared.ko_colors
nk = len(ks)
beats = list(range(n_beats))
data = []
log_vars = [time, vm, nai, inak, inaca, ina, inal, inab]

for k in ks:
    print(f'Simulating {k}mM')

    if debug:
        k = 2.5

    ts = []
    aps = []
    nais = []
    qinas = []
    qinabs = []
    qinals = []
    qinaks = []
    qinacas = []

    s.reset()
    s.set_tolerance(1e-8, 1e-8)
    s.set_constant(ko, k)
    for beat in beats:
        if beat % 100 == 0:
            print('.', end='')
        d = s.run(1000, log=log_vars).npview()
        ts.append(d[time])
        aps.append(d[vm])
        nais.append(d[nai])
        # pA * 1s = pC, *1e-3 = nC
        qinas.append(d.integrate(ina) * 1e-3)
        qinabs.append(d.integrate(inab) * 1e-3)
        qinals.append(d.integrate(inal) * 1e-3)
        qinaks.append(d.integrate(inak) * 1e-3)
        qinacas.append(d.integrate(inaca) * 1e-3)
    print('')

    data.append([ts, aps, nais, qinas, qinabs, qinals, qinaks, qinacas])
    del ts, aps, nais, qinas, qinabs, qinals, qinaks, qinacas

    if debug:
        for k in range(len(ks) - 1):
            data.append(data[0])
        break

#
# Create figure
#
fig = plt.figure(figsize=(9, 8.5))  # Two-column size
fig.subplots_adjust(0.075, 0.06, 0.985, 0.91)
grid = matplotlib.gridspec.GridSpec(4, 4, wspace=0.60, hspace=0.50)

# Add model name
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.980, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.003, 0.895, 'A', letter_font)
fig.text(0.505, 0.895, 'B', letter_font)
fig.text(0.003, 0.670, 'C', letter_font)
fig.text(0.510, 0.670, 'D', letter_font)
fig.text(0.003, 0.440, 'E', letter_font)

#
# APs
#
aplim = -80, 30
ax1 = fig.add_subplot(grid[0, 0])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('V (mV)')
ax1.set_ylim(aplim)
ax2 = fig.add_subplot(grid[0, 1])
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('V (mV)')
ax2.set_ylim(aplim)
for k, c, d in zip(ks, cs, data):
    ts, aps = d[0], d[1]
    ax1.plot(ts[1] * 1e-3, aps[1], label=f'{k} mM', color=c)
    ax2.plot(ts[-1] * 1e-3, aps[-1], label=f'{k} mM', color=c)
ax1.text(0.60, 0.90, text_1, transform=ax1.transAxes)
ax2.text(0.60, 0.90, text_2, transform=ax2.transAxes)
ax1.legend(loc=(1, 1.1), ncol=4)

#
# [Na]i
#
ax6 = fig.add_subplot(grid[0, 2:])
ax6.set_xlabel('Time (s)')
ax6.set_ylabel('$[Na^+]_i$ (mM)')
for k, c, d in zip(ks, cs, data):
    nas = d[2]
    na = np.array([x[0] for x in nas])
    ax6.plot(beats, na, label=f'{k} mM', color=c, zorder=99)

axi = ax6.inset_axes([0.06, 0.23, 0.25, 0.4])
axi.tick_params(labelsize='small')
axi.set_xticklabels([])
axi.set_yticklabels([])
axi.set_ylim(-60, 10)
i = 90
axi.set_xlim(i * cl, i * cl + 500)
axi.plot(data[0][0][i], data[0][1][i], color='tab:orange')
ax6.arrow(i, 9.09, -40, -0.2, zorder=99, color='tab:orange')

axi = ax6.inset_axes([0.54, 0.23, 0.25, 0.4])
axi.tick_params(labelsize='small')
axi.set_xticklabels([])
axi.set_yticklabels([])
axi.set_ylim(-60, 10)
i = 180
axi.set_xlim(i * cl, i * cl + 500)
axi.plot(data[0][0][i], data[0][1][i], color='tab:purple')
ax6.arrow(i + 6, 8.77, 54, 0.08, zorder=99, color='tab:purple')

#
# Na availability
#
ax3 = fig.add_subplot(grid[1, 0])
ax3.set_ylabel('Available $I_{Na}$')
ax3.set_xlabel('V (mV)')
ax3.set_xlim(-95, -40)
ax3.set_ylim(-0.05, 1.05)
ax3.plot(vs, ss, color='#999999')
ax4 = fig.add_subplot(grid[1, 1])
ax4.set_ylabel('Available $I_{Na}$')
ax4.set_xlabel('V (mV)')
ax4.set_xlim(-95, -40)
ax4.set_ylim(-0.05, 1.05)
ax4.plot(vs, ss, color='#999999')
for k, c, d in zip(ks, cs, data):
    aps = d[1]
    vr1 = aps[1][0]
    vr2 = aps[-1][0]
    ax3.plot([vr1], [f(vr1)], 'o', markersize=8, fillstyle='none', color=c)
    ax4.plot([vr2], [f(vr2)], 'o', markersize=8, fillstyle='none', color=c)
ax3.text(0.60, 0.90, text_1, transform=ax3.transAxes)
ax4.text(0.60, 0.90, text_2, transform=ax4.transAxes)

#
# Vr
#
ax5 = fig.add_subplot(grid[1, 2:])
ax5.set_xlabel('Time (s)')
ax5.set_ylabel('$V_r$ (mV)')
for k, c, d in zip(ks, cs, data):
    aps = d[1]
    vr = np.array([x[0] for x in aps])
    ax5.plot(beats, vr, label=f'{k} mM', color=c)


#
# Na ions carried per beat, over time
#
axl = []
qlim = 38, 93
sg = matplotlib.gridspec.GridSpecFromSubplotSpec(
    1, 4, subplot_spec=grid[2:, 0:4], hspace=0, wspace=0.05)
for i, k in enumerate(ks):
    ts, aps, nais, qinas, qinabs, qinals, qinaks, qinacas = data[i]

    ax = fig.add_subplot(sg[0, len(ks) - 1 - i])
    axl.append(ax)
    ax.text(n_beats // 2, 91, f'{k} mM', horizontalalignment='center')

    ax.set_xlabel('Time (s)')
    if i == 3:
        ax.set_ylabel('$Na^+$ carried per cycle (nC)')
    else:
        ax.set_yticklabels([])
    ax.set_ylim(*qlim)
    qinaks1 = 3 * np.array([x[-1] for x in qinaks])
    qinacas1 = -3 * np.array([x[-1] for x in qinacas])
    qinas1 = -1 * np.array([x[-1] for x in qinas])
    qinabs1 = -1 * np.array([x[-1] for x in qinabs])
    qinals1 = -1 * np.array([x[-1] for x in qinals])

    cmap = matplotlib.cm.get_cmap('tab10')
    ax.plot(beats, qinaks1, label='$I_{NaK}$', color=cmap(0))

    out = qinacas1
    ax.fill_between(beats, out * 0, out, color=cmap(1), alpha=0.1)
    ax.plot(beats, out, label='$I_{NaCa}$', color=cmap(1))

    plus = out + qinabs1
    ax.fill_between(beats, out, plus, color=cmap(2), alpha=0.1)
    out = plus
    ax.plot(beats, out, label='... + $I_{NaB}$', color=cmap(2))

    plus = out + qinas1
    ax.fill_between(beats, out, plus, color=cmap(3), alpha=0.1)
    out = plus
    ax.plot(beats, out, label='... + $I_{Na}$', color=cmap(3))

    if i == 3:
        ax.legend()

fig.align_ylabels([ax1, ax5, axl[-1]])

# Show / store
path = os.path.join(fpath, fname)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
plt.savefig(path + '.jpg', dpi=600)
print('Done')
