#!/usr/bin/env python3
#
# Comparison of Grandi model currents to data by Van Wagoner et al.
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
fname = 'figure-d1-restitution'
debug = 'debug' in sys.argv

# Steady-state files
fdata = 'states'

# Simulate?
simulate = (len(sys.argv) > 1)

if simulate:

    # Load model
    m = shared.model()
    name = m.name()

    # Get some variables
    shared.prepare_model(m, pre_pace=False)
    v = m.labelx('membrane_potential')
    g = m.labelx('I_CaL_scaling')

    # Run simulations
    max_beats = 300000
    rel_tol = 1e-4
    cls = [500, 600, 750, 1000, 2000]
    s = myokit.Simulation(m)
    s.set_tolerance(1e-8, 1e-8)
    d1s = []
    d2s = []
    for cl in cls:
        p = myokit.pacing.blocktrain(cl, duration=0.5, offset=0)

        print(f'Pre-pacing model 1 at cl={cl}ms')
        scale = 1
        s.set_constant(g, scale)
        if debug:
            x = m.state()
        else:
            m.get(g).set_rhs(scale)
            path = os.path.join(fdata, f'{name}-cl-{cl}-of-0-ical-{scale}.txt')
            if 'koivumaki' in name and os.path.isfile(path):
                x = myokit.load_state(path)
            else:
                x = shared.limit_cycle(
                    m, p, max_beats=max_beats, rel_tol=rel_tol, path=path)
        s.reset()
        s.set_protocol(p)
        s.set_state(x)
        d1s.append(s.run(min(1000, cl)))

        if debug:
            for i in range(4):
                d1s.append(d1s[0])

        print(f'Pre-pacing model 1 at cl={cl}ms with half ICaL')
        scale = 0.5
        s.set_constant(g, scale)
        if debug:
            x = m.state()
        else:
            m.get(g).set_rhs(scale)
            path = os.path.join(fdata, f'{name}-cl-{cl}-of-0-ical-{scale}.txt')
            if 'koivumaki' in name and os.path.isfile(path):
                x = myokit.load_state(path)
            else:
                x = shared.limit_cycle(
                    m, p, max_beats=max_beats, rel_tol=rel_tol, path=path)
        s.reset()
        s.set_protocol(p)
        s.set_state(x)
        d2s.append(s.run(cl))

        if debug:
            for i in range(4):
                d2s.append(d2s[0])
            break

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
    d['x'] -= d['x'][0]
    d['y'] -= zero_a if x == 'a' else zero_b
    return d

# Create figure
cols = 2 if simulate else 1
fig = plt.figure(figsize=(9, 6)) # Two-column size
fig.subplots_adjust(0.07, 0.08, 0.98, 0.95)
grid = matplotlib.gridspec.GridSpec(3, cols, hspace=0.5, wspace=0.25)

# Add panel letters
letter_font = {'fontsize': 12, 'horizontalalignment': 'center'}
fig.text(0.265, 0.97, 'Van Wagoner et al. 1999', letter_font)
if simulate:
    fig.text(0.770, 0.97, shared.fancy_name(m), letter_font)

# AP, Control
d500 = csv('a', 500)
d600 = csv('a', 600)
d750 = csv('a', 750)
d1000 = csv('a', 1000)
d2000 = csv('a', 2000)

ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel('Time since upstroke (ms)')
ax.set_ylabel('V (mV)')
ax.set_xlim(-50, 610)
ax.set_ylim(-89, 41)
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40])
ax.axhline(0, ls=':', color='#aaaaaa')
ax.plot(d500['x'], d500['y'], label='500ms, Ctrl')
ax.plot(d600['x'], d600['y'], label='600ms')
ax.plot(d750['x'], d750['y'], label='750ms')
ax.plot(d1000['x'], d1000['y'], label='1000ms')
ax.plot(d2000['x'], d2000['y'], label='2000ms')
ax.legend()

if simulate:
    ax = fig.add_subplot(grid[0, 1])
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('V (mV)')
    ax.set_xlim(-50, 610)
    ax.set_ylim(-89, 41)
    ax.set_yticks([-80, -60, -40, -20, 0, 20, 40])
    ax.axhline(0, ls=':', color='#aaaaaa')
    ax.plot(d1s[0].time(), d1s[0][v])
    ax.plot(d1s[1].time(), d1s[1][v])
    ax.plot(d1s[2].time(), d1s[2][v])
    ax.plot(d1s[3].time(), d1s[3][v])
    ax.plot(d1s[4].time(), d1s[4][v])

# AP, Nifedepine
d500 = csv('b', 500)
d600 = csv('b', 600)
d750 = csv('b', 750)
d1000 = csv('b', 1000)
d2000 = csv('b', 2000)

ax = fig.add_subplot(grid[1, 0])
ax.set_xlabel('Time since upstroke (ms)')
ax.set_ylabel('V (mV)')
ax.set_xlim(-20, 610)
ax.set_ylim(-89, 41)
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40])
ax.axhline(0, ls=':', color='#aaaaaa')
ax.plot(d500['x'], d500['y'], label='500ms, Nif')
ax.plot(d600['x'], d600['y'], label='600ms')
ax.plot(d750['x'], d750['y'], label='750ms')
ax.plot(d1000['x'], d1000['y'], label='1000ms')
ax.plot(d2000['x'], d2000['y'], label='2000ms')
ax.legend()

if simulate:
    ax = fig.add_subplot(grid[1, 1])
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('V (mV)')
    ax.set_xlim(-20, 610)
    ax.set_ylim(-89, 41)
    ax.set_yticks([-80, -60, -40, -20, 0, 20, 40])
    ax.axhline(0, ls=':', color='#aaaaaa')
    ax.plot(d2s[0].time(), d2s[0][v])
    ax.plot(d2s[1].time(), d2s[1][v])
    ax.plot(d2s[2].time(), d2s[2][v])
    ax.plot(d2s[3].time(), d2s[3][v])
    ax.plot(d2s[4].time(), d2s[4][v])
    ax.text(400, 20, '50% ICaL block')

# APD
d1 = myokit.DataLog.load_csv(
    os.path.join(dpath, 'van-wagoner-figure-5c-APD90Ctrl.csv'))
d2 = myokit.DataLog.load_csv(
    os.path.join(dpath, 'van-wagoner-figure-5c-APD90Nif.csv'))
d3 = myokit.DataLog.load_csv(
    os.path.join(dpath, 'van-wagoner-figure-5c-APD50Ctrl.csv'))
d4 = myokit.DataLog.load_csv(
    os.path.join(dpath, 'van-wagoner-figure-5c-APD50Nif.csv'))

ax = fig.add_subplot(grid[2, 0])
ax.set_xlabel('Cycle length (ms)')
ax.set_ylabel('Duration (ms)')
ax.set_xlim(400, 2100)
ax.set_ylim(0, 440)
ax.plot(d1['x'], d1['y'], 'ks-', label='APD90, Ctrl')
ax.plot(d2['x'], d2['y'], 'ks-', markerfacecolor='w', label='APD90, Nif')
ax.plot(d3['x'], d3['y'], 'ko-', label='APD50, Ctrl')
ax.plot(d4['x'], d4['y'], 'ko-', markerfacecolor='w', label='APD50, Nif')
ax.legend(loc=(0.5, 0.2)).get_frame().set_alpha(1)

if simulate:
    ax = fig.add_subplot(grid[2, 1])
    ax.set_xlabel('Cycle length (ms)')
    ax.set_ylabel('Duration (ms)')
    ax.set_xlim(400, 2100)
    ax.set_ylim(0, 440)
    ax.plot(cls, apd90_c, 'ks-', label='APD90, Ctrl')
    ax.plot(cls, apd90_n, 'ks-', markerfacecolor='w', label='APD90, Nif')
    ax.plot(cls, apd50_c, 'ko-', label='APD50, Ctrl')
    ax.plot(cls, apd50_n, 'ko-', markerfacecolor='w', label='APD50, Nif')
    #ax.legend(loc=(0.55, 0.2)).get_frame().set_alpha(1)


# Show / store
if not simulate:
    name = 'data-only'
path = os.path.join(fpath, fname + '-' + name)
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
