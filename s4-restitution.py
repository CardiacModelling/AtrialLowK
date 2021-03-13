#!/usr/bin/env python3
#
# Comparison of resitutation characteristics to data from Van Wagoner et al.
#
# Van Wagoner, D. R., Pond, A. L., Lamorgese, M., Rossie, S. S., McCarthy,
# P. M., & Nerbonne, J. M. (1999). Atrial L-type Ca2+ currents and human atrial
# fibrillation. Circulation research, 85(5), 428-436.
#
import matplotlib.pyplot as plt
import matplotlib.gridspec
import myokit
import numpy as np
import os
import shared
import sys

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-s4-restitution'
debug = 'debug' in sys.argv

# Simulate?
simulate = (len(sys.argv) > 1)

cls = [500, 600, 750, 1000, 2000]
max_beats = 30000
if simulate:

    # Load model
    m = shared.model()
    shared.splash(fname, m)
    name = m.name()

    # Prepare model, leaving cleft space alone
    shared.prepare_model(m, fix_cleft_ko=False, pre_pace=False)
    rel_tol = 0.6 if 'koiv' in m.name() else 1e-5

    # Get some variables
    v = m.labelx('membrane_potential')
    g = m.labelx('I_CaL_scaling')

    # Run simulations
    s = myokit.Simulation(m)
    s.set_tolerance(1e-8, 1e-8)
    d1s = []
    d2s = []
    for cl in cls:
        # Create protocol, offset 0 to match experiments
        p = myokit.pacing.blocktrain(cl, duration=0.5, offset=0)

        print(f'Pre-pacing model 1 at cl={cl}ms')
        scale = 1
        m.get(g).set_rhs(scale)
        x, _ = shared.cached_limit_cycle(
            m, p, rel_tol=rel_tol, max_beats=max_beats)

        s.reset()
        s.set_constant(g, scale)
        s.set_protocol(p)
        s.set_state(x)
        d1s.append(s.run(min(1000, cl)))

        print(f'Pre-pacing model 1 at cl={cl}ms with half ICaL')
        scale = 0.5
        m.get(g).set_rhs(scale)
        x, _ = shared.cached_limit_cycle(
            m, p, rel_tol=rel_tol, max_beats=max_beats,
            suffix=f'-ical-{shared.fstr(scale)}')

        s.reset()
        s.set_constant(g, scale)
        s.set_protocol(p)
        s.set_state(x)
        d2s.append(s.run(cl))

    # Calculate APDs
    apd50_c = []
    apd50_n = []
    apd90_c = []
    apd90_n = []
    for d1, d2 in zip(d1s, d2s):
        v0 = d1[v][-1]
        apd50_c.append(d1.apd(v.qname(), threshold=0.5 * v0)['duration'][0])
        apd90_c.append(d1.apd(v.qname(), threshold=0.9 * v0)['duration'][0])
        v0 = d2[v][-1]
        apd50_n.append(d2.apd(v.qname(), threshold=0.5 * v0)['duration'][0])
        apd90_n.append(d2.apd(v.qname(), threshold=0.9 * v0)['duration'][0])

else:
    shared.splash(fname)
    name = 'Van Wagoner et al. 1999'

    # Load and adjust literature data
    dpath = 'literature-data'
    zero_a = os.path.join(dpath, 'van-wagoner-figure-5a-zero.csv')
    zero_a = myokit.DataLog.load_csv(zero_a).npview()
    zero_a = np.mean(zero_a['y'])
    zero_b = os.path.join(dpath, 'van-wagoner-figure-5b-zero.csv')
    zero_b = myokit.DataLog.load_csv(zero_b).npview()
    zero_b = np.mean(zero_b['y'])

    def csv(x, y):
        path = os.path.join(dpath, f'van-wagoner-figure-5{x}-{y}.csv')
        d = myokit.DataLog.load_csv(path).npview()
        d.set_time_key('x')
        d['x'] -= d['x'][0]
        d['y'] -= zero_a if x == 'a' else zero_b
        return d

    # Voltage key
    v = 'y'

    # AP, Control
    d1s = []
    d1s.append(csv('a', 500))
    d1s.append(csv('a', 600))
    d1s.append(csv('a', 750))
    d1s.append(csv('a', 1000))
    d1s.append(csv('a', 2000))

    # AP, Nifedipine
    d2s = []
    d2s.append(csv('b', 500))
    d2s.append(csv('b', 600))
    d2s.append(csv('b', 750))
    d2s.append(csv('b', 1000))
    d2s.append(csv('b', 2000))

    # APD
    apd50_c = myokit.DataLog.load_csv(
        os.path.join(dpath, 'van-wagoner-figure-5c-APD90Ctrl.csv'))['y']
    apd50_n = myokit.DataLog.load_csv(
        os.path.join(dpath, 'van-wagoner-figure-5c-APD90Nif.csv'))['y']
    apd90_c = myokit.DataLog.load_csv(
        os.path.join(dpath, 'van-wagoner-figure-5c-APD50Ctrl.csv'))['y']
    apd90_n = myokit.DataLog.load_csv(
        os.path.join(dpath, 'van-wagoner-figure-5c-APD50Nif.csv'))['y']


# Create figure
fig = plt.figure(figsize=(9, 4))  # Two-column size
fig.subplots_adjust(0.07, 0.11, 0.95, 0.98)
grid = matplotlib.gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.25)

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

ax = fig.add_subplot(grid[:, 0])
ax.set_xlabel('Time since upstroke (ms)')
ax.set_ylabel('V (mV)')
ax.set_xlim(-50, 610)
ax.set_ylim(-89, 41)
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40])
ax.axhline(0, ls=':', color='#aaaaaa')
ax.plot(d1s[0].time(), d1s[0][v], label='500ms, Ctrl')
ax.plot(d1s[1].time(), d1s[1][v], label='600ms')
ax.plot(d1s[2].time(), d1s[2][v], label='750ms')
ax.plot(d1s[3].time(), d1s[3][v], label='1000ms')
ax.plot(d1s[4].time(), d1s[4][v], label='2000ms')
ax.legend(loc='upper right')

# AP, Nifedepine

ax = fig.add_subplot(grid[0, 1])
ax.set_xlabel('Time since upstroke (ms)')
ax.set_ylabel('V (mV)')
ax.set_xlim(-20, 610)
ax.set_ylim(-89, 41)
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40])
ax.axhline(0, ls=':', color='#aaaaaa')
ax.plot(d2s[0].time(), d2s[0][v], label='500ms, Nif')
ax.plot(d1s[1].time(), d1s[1][v], label='600ms')
ax.plot(d1s[2].time(), d1s[2][v], label='750ms')
ax.plot(d1s[3].time(), d1s[3][v], label='1000ms')
ax.plot(d1s[4].time(), d1s[4][v], label='2000ms')
ax.text(400, 20, '50% ICaL block' if simulate else 'Nifedipine')

ax = fig.add_subplot(grid[1, 1])
ax.set_xlabel('Cycle length (ms)')
ax.set_ylabel('Duration (ms)')
ax.set_xlim(400, 2100)
ax.set_ylim(0, 440)
ax.plot(cls, apd90_c, 'ks-', label='APD90, Ctrl')
ax.plot(cls, apd90_n, 'ks-', markerfacecolor='w', label='APD90, Nif')
ax.plot(cls, apd50_c, 'ko-', label='APD50, Ctrl')
ax.plot(cls, apd50_n, 'ko-', markerfacecolor='w', label='APD50, Nif')
ax.legend(loc=(0.5, 0.2)).get_frame().set_alpha(1)


# Show / store
if not simulate:
    name = 'data'
path = os.path.join(fpath, fname + '-' + name)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
