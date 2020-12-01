#!/usr/bin/env python3
#
# Model figure 5: Strength-duration curves
#
import os
import sys
import numpy as np
import matplotlib.patches
import matplotlib.pyplot as plt
import myokit
import shared

# Get path for figure
fpath = shared.figure_dir()
fname = 'figure-m5-sd'
debug = 'debug' in sys.argv
printing = False

# Load model
model = shared.model()
name = model.name()

# Steady-state files
fdata = 'states'

# Strength-duration files
sdata = 'sd-data'
if not os.path.isdir(sdata):
    os.makedirs(sdata)

# Cycle length
cl = 1000

# Time to pre-pace after changing [K]o level
n_beats = 1
if 'long' in sys.argv:
    n_beats = 15 * 60
print(f'Pacing {n_beats} beats after Ko change.')

# Prepare model
shared.prepare_model(model, fix_cleft_ko=True, pre_pace=False)

# Create protocol (with stimulus at t=0) for pre-pacing
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=0)

# Potassium levels
ks = shared.ko_levels
cs = shared.ko_colors
if debug:
    ks, cs = ks[2:3], cs[2:3]

# Get some variables
tt = model.time()
vm = model.labelx('membrane_potential')
st = model.labelx('stimulus_current')
ko = model.labelx('K_o')

# Pre-pace
print(f'Pre-pacing model at cl={cl}ms')
if debug:
    x = model.state()
else:
    # Note: cfix indicates cleft-fix, differs from restitution figure where
    # cleft concentration can vary!
    path = os.path.join(fdata, f'{name}-cl-{cl}-of-0-ical-1-cfix.txt')
    if 'koivumaki' in name and os.path.isfile(path):
        x = myokit.load_state(path)
    else:
        max_beats = 300000
        rel_tol = 1e-4
        x = shared.limit_cycle(
            model, protocol, max_beats=max_beats, rel_tol=rel_tol, path=path)
model.set_state(x)

# Get typical stimulus amplitude
paced = model.clone()
x = paced.binding('pace')
x.set_binding(None)
x.set_rhs(1)
default_amplitude = -paced.labelx('stimulus_current').eval()
print(f'Default amplitude: -{default_amplitude} A/F')
del(paced, x)

# Update i_stim formulation
sa = st.parent().add_variable_allow_renaming('amp')
sa.set_unit('A/F')
sa.set_rhs(default_amplitude)
st.set_rhs(myokit.Multiply(
    myokit.PrefixMinus(myokit.Name(sa)), myokit.Name(model.bindingx('pace'))))

# Durations to test
#durations = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
durations = np.logspace(-2, 1, 100)
ds = []

# Start testing
for k in ks:
    print(f'Testing {k} mM')
    k_str = str(float(k)).replace('.', '-')
    path = os.path.join(sdata, f'sd-{name}-ko-{k_str}-b-{n_beats}.csv')
    try:
        ds.append(myokit.DataLog.load_csv(path))
        print(f'Loaded data from {path}')
        continue
    except Exception:
        pass

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    s.set_constant(ko, k)

    # Pre-pace for 1s
    if not debug:
        print(f'Pre-pacing for {n_beats * cl} ms')
        s.pre(n_beats * cl)

    # Get an expected Vmax and APD
    d = s.run(cl, log=[tt, vm]).npview()
    vmin = np.min(d[vm])
    vmax = np.max(d[vm])
    vacc = -30  # Minimum accepted max V
    vap = 0.9 * vmin

    try:
        apd = d.apd(threshold=vap)['duration'][0]
    except IndexError:
        apd = float('nan')
    tmax = min(2 * apd, cl)
    print(f'vmin = {vmin}')
    print(f'vmax = {vmax}')
    print(f'APD = {apd}')

    # Method to test if an AP occurred
    def test(s, a):
        if printing:
            print(f'  Testing {a} A/F')
        s.reset()
        s.set_constant(sa, a)

        # Check if it depolarises
        d = s.run(10, log=[tt, vm])
        if np.max(d[vm]) < vacc:
            if printing:
                print('   No depolarisation')
            return False

        # Check if we get an acceptable APD
        d = s.run(tmax, log=d)
        a = d.apd(threshold=vap)['duration']
        if len(a) != 1:
            if printing:
                print('   No APD')
            return False
        if a[0] < apd * 0.5:
            if printing:
                print('   APD too short')
            return False
        return True

    # Test each duration
    amplitudes = []
    if np.isnan(apd):
        amplitudes = [apd] * len(durations)
    else:
        for duration in durations:
            print(f' Testing {duration} ms')

            # Create protocol
            s.set_protocol(myokit.pacing.blocktrain(cl, duration))

            # Search for minimum stimulus
            amin = 0
            amax = default_amplitude / 8

            amax_found = False
            for i in range(10):
                try:
                    if test(s, amax):
                        amax_found = True
                        break
                except myokit.SimulationError:
                    pass
                amax *= 2
            if not amax_found:
                print('  Failed to find upper bound for search')
                amplitudes.append(float('nan'))
                continue

            # Zero must lie in between. Start bisection search
            failures = 0
            while abs(amax - amin) > 1e-2:
                a = 0.5 * (amax + amin)
                try:
                    if test(s, a):
                        amax = a
                    else:
                        amin = a
                except myokit.SimulationError:
                    print('  Simulation failed!')
                    failures += 1
                    if failures == 10:
                        a = float('nan')
                        break
                    else:
                        amax += 0.01 * (amax - amin)

            amplitudes.append(a)

    # Store
    d = myokit.DataLog()
    d['duration'] = durations
    d['amplitude'] = amplitudes
    d.save_csv(path)
    ds.append(d)
print('Done')


#
# Create figure
#
fig = plt.figure(figsize=(9, 5))  # Two-column size
fig.subplots_adjust(0.08, 0.12, 0.95, 0.98)

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

# Plot curves
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('Stimulus duration (ms)')
ax.set_ylabel('-1x Stimulus amplitude (A/F)')
ax.axvline(0.5, color='#cccccc')

for k, c, d in zip(ks, cs, ds):
    ax.plot(d['duration'], d['amplitude'], color=c)
ax.set_xlim(0, 10)
ax.set_ylim(0, 100)

if n_beats == 1:
    atext = 'After 1 second'
else:
    atext = n_beats / 60
    atext = int(atext) if int(atext) == atext else round(atext, 1)
    atext = f'After {atext} minutes'
ax.text(0.98, 0.95, atext, transform=ax.transAxes, ha='right')

# Show / store
path = os.path.join(fpath, f'{fname}-{name}-b-{n_beats}')
print('Saving to ' + path)
plt.savefig(path + '.png')
plt.savefig(path + '.pdf')
print('Done')
