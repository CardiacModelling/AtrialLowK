#!/usr/bin/env python3
#
# Sets up strand simulations and
#  - Pre-paces for a number of beats
#  - Applies a Ko change and simulates N extra beat(s).
#  - Simulates a beat to calculate CV.
#  - Rewinds 1 beat, then finds the refractory period.
#
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import myokit
import shared

# Get path for figures
fpath = shared.figure_dir()

# Test beat to find g and f?
test_beat = False
try:
    if sys.argv[2] in 'fg':
        test_beat = sys.argv[2]
        print('Running test beat.')
except IndexError:
    pass

# Parse command line arguments
level = 5.4
cvt = 5
af = False
long_term = False
long_text = ''
if not test_beat:
    # Get ko level and CV threshold for a "propagating wave"
    try:
        level = float(sys.argv[2])
        cvt = float(sys.argv[3])
    except (IndexError, ValueError):
        print('Syntax: ' + sys.argv[0] + ' model ko cvt (af)')
        sys.exit(1)
    cvt = int(cvt) if int(cvt) == float(cvt) else cvt

    # Get AF mode
    try:
        af = (sys.argv[4].lower() == 'af')
    except IndexError:
        pass

    # Get long run
    if len(sys.argv) > 4 and sys.argv[-1] == 'long':
        long_term = True
        long_text = '-long'

# Get path for data
dpath = 'strand-data'

# Load model
model = shared.model()
name1 = model.name() + ('-af' if af else '')
name2 = f'{name1}-cvt-{cvt}-ko-' + str(level).replace('.', '-')

# Show details
print('=' * 79)
print(f'Running for {model.name()} with [K]o-change to {level}')
print(f'Propagating wave threshold: CV > {cvt}')
print(f'Running {"LOOOOONG 15min" if long_term else "1s"} simulation')

# Set af model
if af:
    print('Enabling Chronic AF mode')
    shared.enable_af(model)
else:
    print('Chronic AF mode not enabled')

# Stimulus protocol timing
cl = 1000
t0 = 50
sd = 5
print(f'Stimulus of {sd} ms applied every {cl} ms with a {t0} offset.')

# Number of pre-pacing beats
# WARNING: Not implemented everywhere
n_pre = 50

#
# Note: The code below is almost independent of the choice of variables above.
# However: It assumes that n_pre is constant, and refers to simulations with
# some number of n_pre and some number of n_beats as sim_{n_pre+n_beats}. But
# if n_pre varies, this becomes blurry / unidentifiable.
# TODO: Either re-write as having n_pre fixed, or change the file format to
# differentiate between n_pre and n_beats.
#

# Prepare model
shared.prepare_model(model, fix_cleft_ko=True)

# Get some variables
tt = model.time()
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')

# Ensure suitable single-cell state is set (hardcoded for baseline)
use_cached = True
if af:
    # Pre-pace to periodic orbit with AF mode enabled.
    fname = name1 + '-single-cell-state.txt'
    path = os.path.join(dpath, fname)
    use_cached &= os.path.isfile(path)
    if use_cached:
        print('Loading pre-paced single-cell AF state')
        model.set_state(myokit.load_state(path))
    else:
        print(f'Pre-pacing: {name1}')
        p = myokit.pacing.blocktrain(cl, duration=0.5, offset=t0)
        model.set_state(shared.limit_cycle(model, p, max_beats=50000))
        print(model.format_state(model.state()))
        del(p)
        print('Storing single-cell AF state')
        myokit.save_state(path, model.state())

# Conduction velocity
cvargs = dict(name=vm.qname(), length=0.01)   # Cell length in cm

# Set conductance and stimulus level
# (pA/pF)/mV = nS/pF = S/mF
# Note: Conductance was set first, then minimum stimulus, but the two are
# slightly dependent.
if model.name() == 'courtemanche-1998':
    g = 26.98
    f = 0.171
elif model.name() == 'modified-grandi-2011':
    g = 23.32
    f = 0.168
elif model.name() == 'grandi-2011':
    g = 30.47
    f = 0.237
elif model.name() == 'koivumaki-2011':
    g = 36.09
    f = 0.198
elif model.name() == 'maleckar-2008':
    g = 38.16
    f = 0.198
elif model.name() == 'ni-2017':
    g = 23.20
    f = 0.170
elif model.name() == 'nygren-1998':
    g = 41.95
    f = 0.215
elif model.name() == 'voigt-heijman-2013':
    g = 11.29
    f = 0.140
else:
    raise NotImplementedError(f'No conductance known for {name1}.')

# Set stimulus amplitude as multiple of minimum
if test_beat != 'f':
    factor = 2
    print(f'Using stimulus that is {factor} times the minimum amplitude.')
    f *= factor
    del(factor)

# Create protocol
pargs = dict(period=cl, level=f, duration=sd)
protocol = myokit.pacing.blocktrain(offset=t0, **pargs)

# Create simulation
s = myokit.SimulationOpenCL(
    model, protocol, ncells=(200), rl=True, precision=myokit.DOUBLE_PRECISION)
s.set_step_size(0.001)
s.set_conductance(g)

# Test stimulus size and/or conductance
if test_beat:
    print('Simulating test beat for '
          + ('conductance' if test_beat == 'g' else 'stimulus amplitude'))
    print(f'f={f}, g={g}')
    d = s.run(100)
    print('Done')
    b = d.block1d()
    cv = b.cv(**cvargs)
    print(f'CV: {cv} cm/s')
    sys.exit(1)

# Pre-pace for n_pre beats at baseline ko
fname = f'{name1}-state-{n_pre}.zip'
path = os.path.join(dpath, fname)
use_cached &= os.path.isfile(path)
if use_cached:
    print(f'Loading state after {n_pre} beats at baseline Ko')
    state_pre = myokit.load_state_bin(path)
    s.set_state(state_pre)
else:
    print(f'Pre-pacing {n_pre} beats at baseline Ko')
    print(f'0        1         2         3         4         5')
    print(f'12345678901234567890123456789012345678901234567890')
    for i in range(n_pre):
        s.run(cl, log=myokit.LOG_NONE)
        if (1 + i) % 50 == 0:
            print(f'. ({1 + i})')
        else:
            print('.', end='')
            sys.stdout.flush()
    print('Storing state')
    state_pre = s.state()
    myokit.save_state_bin(path, state_pre)
    s.set_time(0)
del(fname, path)
print(f'Now at state after {n_pre} baseline beats.')
print(f'In-simulation time: t={s.time()} ms')

# Change [K]o
print(f'Changing [K]o to {level}')
s.set_constant(ko, level)

# Simulate 1st beat at new Ko
n_beats = 1                 # 1 second, assuming CL=1000ms
i_ready = n_pre + n_beats   # Index of pre-paced beat
i_final = i_ready + 1       # Index of test beat
fname1 = f'{name2}-state{i_ready}.zip'
fname2 = f'{name2}-beat{i_ready}.bin'
path1 = os.path.join(dpath, fname1)
path2 = os.path.join(dpath, fname2)
use_cached &= os.path.isfile(path1)
use_cached &= os.path.isfile(path2)
if use_cached:
    print(f'Loading state {i_ready} (state after {n_beats} beat(s) at new Ko)')
    state_ready = myokit.load_state_bin(path1)
    time_ready = n_beats * cl
    s.set_state(state_ready)
    s.set_time(time_ready)
else:
    print(f'Simulating {n_beats} beat(s) pre-pacing at new Ko')
    d = s.run(cl, log=[tt, vm])
    print(f'Storing state {i_ready}')
    state_ready = s.state()
    time_ready = n_beats * cl
    myokit.save_state_bin(path1, state_ready)
    print(f'Storing beat {i_ready}')
    d.block1d().save(path2)
    del(d)
del(fname1, fname2, path1, path2)
print(f'Now at state after {n_pre} baseline beat(s)'
      f' + {n_beats} beat(s) with altered Ko.')
print(f'In-simulation time: t={s.time()} ms')

# Simulate extra beats at lowered K, if needed
if long_term:
    n_done = n_beats
    n_beats = 15 * 60           # 15 minutes, assuming CL=1000ms
    i_ready = n_pre + n_beats   # Index of pre-paced beat
    i_final = i_ready + 1       # Index of test beat
    n_todo = n_beats - n_done

    fname1 = f'{name2}-state{i_ready}.zip'
    path1 = os.path.join(dpath, fname1)
    if use_cached and os.path.isfile(path1):
        print(f'Loading state {i_ready}'
              f' (state after {n_beats} beat(s) at new Ko)')
        state_ready = myokit.load_state_bin(path1)
        time_ready = n_beats * cl
        s.set_state(state_ready)
        s.set_time(time_ready)
    else:
        # Load intermediary beats, if any
        if use_cached:
            i_lo = (n_pre + n_done) // 50
            i_hi = (n_pre + n_todo) // 50 -1
            for i in range(i_hi, i_lo, -1):
                j = 50 * i
                fname2 = f'{name2}-state{j}.zip'
                path2 = os.path.join(dpath, fname2)
                if os.path.isfile(path2):
                    print(f'Loading state {j}')
                    s.set_state(myokit.load_state_bin(path2))
                    j -= n_pre
                    print(f'Now at state after {n_pre} baseline beat(s)'
                          f' + {j} beat(s) at new Ko')
                    s.set_time(j * cl)
                    print(f'In-simulation time: t={j * cl} ms')
                    n_done = j
                    n_todo = n_beats - j
                    break
            del(i_lo, i_hi)
        # At this point, we have loaded an intermediate state, but still need
        # to simulate some more, so use_cached for future files should be
        # false.
        use_cached = False

        # Simulate remaining beats
        print(f'Simulating {n_todo} beats pre-pacing,'
              f' for a total of {n_beats} beats at new Ko.')
        print('12345678901234567890123456789012345678901234567890')
        for i in range(n_todo):
            s.run(cl, log=myokit.LOG_NONE)
            # Save intermediary state
            m = '.'
            j = n_pre + n_done + i + 1
            if (j % 50 == 0) and (i + 1 < n_todo):
                m = 's'
                fname2 = f'{name2}-state{j}.zip'
                path2 = os.path.join(dpath, fname2)
                myokit.save_state_bin(path2, s.state())
            # Show status
            if (1 + i) % 50 == 0:
                print(f'{m} ({1 + i})')
            else:
                print(f'{m}', end='')
                sys.stdout.flush()
        print()
        print(f'Storing state {i_ready}')
        state_ready = s.state()
        time_ready = n_beats * cl
        myokit.save_state_bin(path1, state_ready)

        del(fname2, path2)
    del(fname1, path1, n_done, n_todo)
    print(f'Now at state after {n_pre} baseline beats'
          f' + {n_beats} beat(s) with altered Ko.')
    print(f'In-simulation time: t={s.time()} ms')

# Simulate extra beat to get CV, APD, APD start, and to plot
fname = f'{name2}-beat{i_final}-cv.bin'
path = os.path.join(dpath, fname)
use_cached &= os.path.isfile(path)
if use_cached:
    print(f'Loading beat {i_final} for CV, APD, and AP start')
    b_final_cv = myokit.DataBlock1d.load(path)
    d_final_cv = b_final_cv.to_log()
else:
    print(f'Simulating beat {i_final} for CV, APD, and AP start')
    d_final_cv = s.run(cl, log=[tt, vm])

    print(f'Storing beat {i_final}')
    b_final_cv = d_final_cv.block1d()
    b_final_cv.save(path)

    # Reset simulation to state after Nth beat
    print(f'Resetting to state after beat {i_ready}')
    s.set_time(time_ready)
    s.set_state(state_ready)
del(fname, path)
print(f'Now at state after {n_pre} baseline beats'
      f' + {n_beats} beat(s) with altered Ko.')
print(f'In-simulation time: t={s.time()} ms')

# Calculate CV
cv = b_final_cv.cv(**cvargs)

# Estimate APD
apd_cell = 0
vmin = d_final_cv[vm.qname(), apd_cell][0]
vmax = 40
apd_threshold = vmin + 0.1 * (vmax - vmin)
print(f'Calculating APD in cell {apd_cell}'
      f' with threshold of {apd_threshold} mV')
apd = d_final_cv.apd(f'{apd_cell}.{vm.qname()}', threshold=apd_threshold)
start, apd = apd['start'][0], apd['duration'][0]
del(vmin, vmax, apd_threshold)

print(f'CV: {cv} cm/s')
print(f'AP starts at {start} ms and lasts until t={start + apd} ms')
print(f'APD: {apd} ms')

# Calculate minimum and maximum expected refractory period
# Minimum: Rounded end of APD90, minus some margin
# Note: `start` includes time_ready
tmin = max(time_ready, 10 * int((start + apd) / 10) - 100)
# Maximum: Start time plus a beat, plus the stimulus offset, plus a margin
tmax = time_ready + cl + t0 + 200
print(f'Search range for FRP: {tmin} to {tmax}')
print(f'FRP accepted if CV > {cvt} ms')

# Find refractory time
fname1 = f'{name2}-beat{i_final}-refractory-period.csv'
fname2 = f'{name2}-beat{i_final}-wavelength.csv'
path1 = os.path.join(dpath, fname1)
path2 = os.path.join(dpath, fname2)
use_cached &= os.path.isfile(path1)
use_cached &= os.path.isfile(path2)
if use_cached:
    print('Loading FRP and WL data from disk')
    csv = myokit.DataLog.load_csv(path1).npview()
    ts = csv['ts']
    vs = csv['vs']
    csv = myokit.DataLog.load_csv(path2)
    rp = csv['rp'][0]   # ms
    wl = csv['wl'][0]   # mm
    del(csv)
else:
    # Estimate FRP
    print('Estimating FRP and WL')
    rp = wl = float('nan')
    ts = [tmax]
    vs = [cv]

    # Method to try out a refractory time
    def refrac(s, t, t0, s0):
        # Reset simulation
        s.set_time(t0)
        s.set_state(s0)

        # Remove any protocol
        s.set_protocol(None)

        # Simulate up to just before the stimulus
        cache = None
        t_final = t - 5
        duration = t_final - t0
        if duration > 10:
            print(f'  Simulating from t={t0} ms to t={t_final} ms'
                  f' (duration {duration} ms).')
            s.run(duration, log=myokit.LOG_NONE)
            t0 = t_final
            cache = (t_final, s.state())

        # Set protocol with beat at time t
        # (pargs contains period, duration, and level)
        s.set_protocol(myokit.pacing.blocktrain(offset=t, limit=1, **pargs))

        # Run until a short while after t
        # At 5cm/s, it would take 1.5cm/(5cm/s) = 0.3s for cell 150 to get
        # activated, so 300ms sim should be long enough for even very slow CVs
        # to be measured (cv method ignores first and last 50 cells).
        # At 10cm/s, 150ms should be enough
        t_final = t + 150
        duration = t_final - t0
        print(f'  Simulating from t={t0} ms to t={t_final} ms for a beat at'
              f' t={t} ms (duration {duration} ms).')
        d = s.run(duration, log=[tt, vm])
        b = d.block1d()
        cv = d.block1d().cv(**cvargs)
        cv = max(0, cv)

        return cv, cache

    if cv > cvt:
        # Simulate up to tmin
        # Keep earliest time and state from time guaranteed to be < FRP
        duration = tmin - time_ready
        print(f'  Simulating up to t={tmin} ms (duration {duration} ms)')
        s.run(duration, log=myokit.LOG_NONE)
        time_lower = tmin
        state_lower = s.state()

        # Run bisection search
        t1 = tmin
        t2 = tmax
        ts = []
        vs = []
        while t2 - t1 > 0.1:
            print(f'   Searching [{t1}, {t2}]')

            # Calculate CV at bisection point
            t = t1 + 0.5 * (t2 - t1)
            v, cache = refrac(s, t, time_lower, state_lower)

            # Update the edges, cache simulation state if possible
            if v > cvt:
                print(f'   CV = {v} << down [{t1}, {t2}]')
                t2 = t
            else:
                print(f'   CV = {v} >> up')
                t1 = t

                # Cache earliest state without refraction.
                # Speeds up searches where the FRP is very very long.
                if cache is not None:
                    print(f'   Updating cached state to t={cache[0]} ms')
                    time_lower = cache[0]
                    state_lower = cache[1]

            # Store time pulse was given, and resulting CV
            ts.append(t)
            vs.append(v)

        # Sort times and cvs
        ts, vs = np.array(ts), np.array(vs)
        i = np.argsort(ts)
        ts, vs = ts[i], vs[i]

        # Refractory time and wavelength
        if len(vs) > 0 and max(vs) > cvt:
            # RP = Lowest conducting t, minus t start, minus stimulus offset
            rp = t2 - time_ready - t0
            wl = cv * rp * 0.01 # In mm

    # Store results
    csv = myokit.DataLog(time='ts')
    csv['ts'] = ts
    csv['vs'] = vs
    csv.save_csv(path1)

    csv = myokit.DataLog()
    csv['rp'] = [rp]
    csv['wl'] = [wl]
    csv['cv'] = [cv]
    csv['apd'] = [apd]
    csv.save_csv(path2)
    del(csv)

    # Reset simulation to state after pre-pacing beats
    print(f'Resetting to state after beat {i_ready}')
    s.set_time(time_ready)
    s.set_state(state_ready)
del(fname1, fname2, path1, path2)
print(f'Now at state after {n_pre} baseline beats'
      f' + {n_beats} beat(s) with altered Ko.')
print(f'In-simulation time: t={s.time()} ms')

# Show main results
print(f'APD: {apd} ms')
print(f'CV:  {cv} cm/s')
print(f'RP:  {rp} ms')
print(f'WL:  {wl} mm')

# Simulate final beat with minimum refractory beat
if not np.isfinite(rp):
    print(f'No FRP known: NOT simulating beat {i_final} and refractory beat.')
else:
    fname = f'{name2}-beat{i_final}-refrac.bin'
    path = os.path.join(dpath, fname)
    use_cached &= os.path.isfile(path)
    if use_cached:
        print(f'Loading beat {i_final} and refractory beat')
        b_final_refrac = myokit.DataBlock1d.load(path)
    else:
        print(f'Simulating beat {i_final} and refractory beat')

        # Simulate with normal protocol, until tmin
        duration = tmin - time_ready
        print(f'Simulating from t={time_ready} ms to t={tmin} ms'
              f' (duration {duration} ms)')
        s.set_protocol(protocol)
        d = s.run(duration, log=[tt, vm])

        # Pick final time to simulate to
        t_final = time_ready + 800
        while time_ready + t0 + rp + apd + 100 > t_final:
            t_final += 200

        # Simulate with extra stimulus protocol, from tmin to tmax
        t = time_ready + rp + t0
        s.set_protocol(myokit.pacing.blocktrain(offset=t, limit=1, **pargs))
        duration = t_final - tmin
        print(f'Simulating from t={tmin} ms to t={t_final} ms for a beat at'
              f' t={t} ms (duration {duration} ms).')
        d = s.run(duration, log=d)

        print('Storing final and refractory beat')
        b_final_refrac = d.block1d()
        b_final_refrac.save(path)
        del(d)
    del(fname, path)


#
# Create figures
#
print('Creating figures')

#
# AP and CV versus stimulus time
#
fig = plt.figure(figsize=(4.3, 3))  # One-column size
fig.subplots_adjust(0.14, 0.17, 0.80, 0.96)

# Add model name
fig.patches.append(Rectangle((0.93, 0), 0.07, 1, color='#cccccc', fill=True,
    transform=fig.transFigure, zorder=-1))
name_font = {
    'fontsize': 12,
    'rotation': 'vertical',
    'verticalalignment': 'center',
}
fig.text(0.95, 0.5, name1, name_font)

# Graph AP and CV
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel(f'Vm in cell {1 + apd_cell} (mV)')
ax1.set_ylim(-105, 65)
ax1.plot(d_final_cv.time(), d_final_cv[vm.qname(), apd_cell])
ax2 = ax1.twinx()
ax2.set_ylabel('CV (cm/s)')
ax2.set_ylim(0, 100)
ax2.plot(ts, vs, 'x-', color='tab:orange')

fname = f'{name2}-beat{i_final}-refractory-period'
path = os.path.join(fpath, fname)
plt.savefig(path + '.png')

#
# APs throughout the stimulus
#
if np.isfinite(rp):
    fig = plt.figure(figsize=(9, 2.7))  # Two-column size
    fig.subplots_adjust(0.07, 0.17, 0.95, 0.96)

    # Add model name
    fig.patches.append(Rectangle((0.965, 0), 0.035, 1, color='#cccccc',
        fill=True, transform=fig.transFigure, zorder=-1))
    fig.text(0.975, 0.5, name1, name_font)

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Vm (mV)')
    t = b_final_refrac.time()
    v = b_final_refrac.get1d(vm.qname())
    for j in range(20):
        i = 4 + 10 * j
        ax.plot(t, v[:, i], lw=0.7)

    fname = f'{name2}-beat{i_final}-aps'
    path = os.path.join(fpath, fname)
    plt.savefig(path + '.png')

