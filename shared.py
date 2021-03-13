#!/usr/bin/env python3
#
# Shared code for external K investigations
#
import os
import sys

import matplotlib
import matplotlib.pyplot    # Makes .cm exist!
import myokit
import myokit.lib.hh
import numpy as np


# Update matplotlib styles
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
#matplotlib.rcParams['xtick.minor.visible'] = True
#matplotlib.rcParams['ytick.minor.visible'] = True
matplotlib.rcParams['mathtext.default'] = 'regular'

# Extracellular Potassium concentrations
ko_levels = [2.5, 3.2, 4, 5.4]
#ko_colors = ['k', '#5e3c99', '#b2abd2', '#fdb863', '#e66101']
#ko_colors = ['#1d003b', '#6e4ca9', '#b2abd2', '#fdb863', '#e66101']
ko_colors = ['#000000', '#d7191c', '#fdae61', '#2c7bb6']


# Simulation data
sdata = 'sim-data'


# Models
model_names = {
    'voigt': 'voigt-heijman-2013.mmt',
    'grandi': 'grandi-2011.mmt',
    'courtemanche': 'courtemanche-1998.mmt',
    'ni': 'ni-2017.mmt',
    'nygren': 'nygren-1998.mmt',
    'maleckar': 'maleckar-2008.mmt',
    'koivumaki': 'koivumaki-2011.mmt',
}

# Current colors
cmap = matplotlib.cm.get_cmap('tab20')
current_names = [
    'I_Kur',    # 0
    'I_to',     # 1
    'I_Cl,B',   # 2
    'I_CaL',    # 3
    'I_NaCa',   # 4
    'I_Kr',     # 5
    'I_Ks',     # 6
    'I_K1',     # 7
    'I_NaK',    # 8
    'I_Ca,P',   # 9
    'I_Ca,B',   # 10
    'I_Na,B',   # 11
    'I_KACh',   # 12
    'I_ClCa',   # 13
    'I_Kp',     # 14
    'I_f',      # 15
    'I_Na',     # 16
    'I_NaL',    # 17
]


def cached_limit_cycle(model, protocol, path=None, suffix=None,
                       rel_tol=1e-5, max_beats=4000):
    """
    Pre-paces a model to periodic orbit ("steady state"), or loads a periodic
    state from disk.

    Arguments
    ``model``
        The model to pre-pace
    ``protocol``
        The protocol to pre-pace with
    ``path``
        An optional path (without an extension) used in caching, if set this
        will override the default generated path name.
    ``suffix``
        If a suffix is given, this will be appended to the path.

    Returns a tuple ``(state, cached)`` where ``state`` is the new state, and
    ``cached`` is True if a cached result was returned.
    """
    default_path, cl = _path_cl(model, protocol)
    if path is None:
        path = default_path
    if suffix:
        path += suffix
    path1 = path + '-steady.zip'
    path2 = path + '-stabilising.zip'

    # Load cached and stable result
    try:
        state = model.map_to_state(myokit.load_state_bin(path1))
        print(f'Loaded pre-paced state from {path1}')
        return state, True
    except Exception:
        pass

    # Load previous state from file
    try:
        model = model.clone()
        model.set_state(model.map_to_state(myokit.load_state_bin(path2)))
        print(f'Loaded partially pre-paced state from {str(path2)}')
    except Exception:
        pass

    # Pre-pace
    print('Pre-pacing...')
    state, period = limit_cycle(model, protocol, rel_tol, max_beats)

    # Store stable or unstable state
    if period == 1:
        print(f'Storing pre-paced state to {path1}')
        myokit.save_state_bin(path1, state)
    elif period == 0:
        print(f'Storing partially pre-paced state to {path2}')
        myokit.save_state_bin(path2, state)

    print(model.format_state(model.state()))
    return state, False


def currents(model, colours=False):
    """
    Returns an ordered list of transmembrane current variable names, for a
    current plot.
    """
    name = model.name().lower()
    if 'courtemanche-1998' in name:
        currents = {
            'ikur.i_Kur': 0,
            'ito.i_to': 1,
            'ical.i_Ca_L': 3,
            'inaca.i_NaCa': 4,
            'ikr.i_Kr': 5,
            'iks.i_Ks': 6,
            'ik1.i_K1': 7,
            'inak.i_NaK': 8,
            'ipca.i_PCa': 9,
            'ib.i_B_Ca': 10,
            'ib.i_B_Na': 11,
            'ina.i_Na': 16,
        }
    elif 'grandi-2011' in name:
        currents = {
            'ikur.IKur': 0,
            'ito.Ito': 1,
            'iclb.IClB': 2,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'ipca.IpCa': 9,
            'icab.ICaB': 10,
            'inab.INaB': 11,
            'iclca.IClCa': 13,
            'ikp.IKp': 14,
            'ina.INa': 16,
            'inal.INaL': 17,
        }
    elif 'koivumaki' in name:
        currents = {
            'ikur.IKur': 0,
            'it.It': 1,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'icap.ICaP': 9,
            'icab.ICab': 10,
            'inab.INab': 11,
            'if.If': 15,
            'ina.INa': 16,
        }
    elif 'maleckar-' in name:
        currents = {
            'ikur.i_Kur': 0,
            'it.i_t': 1,
            'ical.i_Ca_L': 3,
            'inaca.i_NaCa': 4,
            'ikr.i_Kr': 5,
            'iks.i_Ks': 6,
            'ik1.i_K1': 7,
            'inak.i_NaK': 8,
            'icap.i_CaP': 9,
            'ib.i_B_Ca': 10,
            'ib.i_B_Na': 11,
            'ikach.i_KACh': 12,
            'ina.i_Na': 16,
        }
    elif 'ni-' in name:
        currents = {
            'ikur.IKur': 0,
            'ito.Ito': 1,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'icap.ICap': 9,
            'ibca.IbCa': 10,
            'ibna.IbNa': 11,
            'ina.INa': 16,
        }
    elif 'nygren' in name:
        currents = {
            'isus.i_sus': 0,
            'it.i_t': 1,
            'ical.iCaL': 3,
            'inaca.i_NaCa': 4,
            'ikr.i_Kr': 5,
            'iks.i_Ks': 6,
            'ik1.i_K1': 7,
            'inak.i_NaK': 8,
            'icap.i_CaP': 9,
            'ib.i_B_Ca': 10,
            'ib.i_B_Na': 11,
            'ina.i_Na': 16,
        }
    elif 'voigt' in name:
        currents = {
            'ikur.IKur': 0,
            'ito.Ito': 1,
            'iclb.IClB': 2,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'ipca.IpCa': 9,
            'icab.ICaB': 10,
            'inab.INaB': 11,
            'ikach.IKACh': 12,
            'iclca.IClCa': 13,
            'ikp.IKp': 14,
            'ina.INa': 16,
            'inal.INaL': 17,
        }
    else:
        # Get all currents
        def rec(parent, currents=set()):
            for var in parent.refs_to():
                if 'tot' in var.name():
                    rec(var, currents)
                else:
                    currents.add(var)
            return currents

        currents = rec(model.labelx('cellular_current'))
        currents = [x.qname() for x in currents]
        currents.sort()
        print('\n'.join(currents))
        print(len(currents))
        raise NotImplementedError('Unknown model: ' + model.name())

    if colours:
        colours = [cmap(x) for x in currents.values()]
        currents = list(currents.keys())
        return currents, colours
    return list(currents.keys())


def default_protocol():
    """ Returns the protocol used in most simulations. """
    return myokit.pacing.blocktrain(1000, duration=0.5, offset=50)


def demote(var):
    """
    Changes a state variable to a non-state variable.

    1. Replaces any references to its derivative with an inlined equation
    2. Sets the variable's RHS to its state value
    3. Demotes the variable.
    """
    if not var.is_state():
        return

    # Replace references to dot(x) by inlined equation
    subst = {var.lhs(): var.rhs()}
    for ref in list(var.refs_by()):
        ref.set_rhs(ref.rhs().clone(subst=subst))

    # Demote
    var.set_rhs(var.state_value())
    var.demote()


def enable_af(model):
    """
    Places the given ``model`` into AF mode.
    """
    name = model.name()
    #if name == 'courtemanche':
    if 'grandi-2011' in name:
        model.get('cell.AF').set_rhs(1)
    #elif 'koivumaki' in name:
    #elif 'maleckar' in name:
    elif 'ni-' in name:
        model.get('af.af').set_rhs(1)
    #elif 'nygren' in name:
    elif 'voigt' in name:
        model.get('cell.AF').set_rhs(1)
    else:
        raise NotImplementedError(f'No AF mode available for {name}.')


def fancy_name(model):
    """Returns a nicer name for a model."""
    name = model.name()
    if name == 'grandi-2011':
        return 'Grandi-Pandit-Voigt et al., 2011'
    if 'courtemanche' in name:
        return 'Courtemanche et al., 1998'
    if 'koivumaki' in name:
        return 'Koivumaki et al., 2011'
    if 'maleckar' in name:
        return 'Maleckar et al., 2008'
    if 'ni' in name:
        return 'Ni et al., 2017'
    if 'nygren' in name:
        return 'Nygren et al., 1998'
    if 'voigt' in name:
        return 'Voigt-Heijman et al., 2013'
    raise NotImplementedError('Unknown model: ' + str(name))


def figure_dir():
    """
    Returns the path to a directory for storing figures, depending on the
    script name.
    """
    path = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    if path[:1] == 'f':
        path = 'figure-' + path[1:]
    else:
        path = 'sup-figure-' + path[1:]
    if not os.path.isdir(path):
        os.makedirs(path)
    return path


def fstr(x):
    """Returns a string of number x, with ``-`` instead of a decimal point."""
    return str(int(x) if float(x) == int(x) else float(x)).replace('.', '-')


def ina_availability(model):
    """
    Creates an INa availability curve for the given model, at different [K]o
    levels.

    Returns a tuple ``(v, a, v1s, a1s, v2s, a2s)`` where ``v`` and ``s`` are
    lists of voltages and availability values. The resting potentials for each
    Ko are given in v1s and v2s, while a1s and a2s provide the corresponding
    availability levels.
    """
    model = model.clone()
    protocol = default_protocol()

    # Get the membrane voltage indice
    vm = model.labelx('membrane_potential').indice()

    # Get steady-state inactivation curve
    try:
        # Courtemanche and grandi families use (h * j) for inactivation
        ina_j_inf = model.labelx('ina_j_inf')
        ina_h_inf = model.labelx('ina_h_inf')
        fj = ina_j_inf.pyfunc()
        fh = ina_h_inf.pyfunc()
        f = lambda v: fj(v) * fh(v)
    except myokit.IncompatibleModelError:
        # Nygren family uses (0.9 * h1 + 0.1 * h2) for inactivation
        ina_h12_inf = model.labelx('ina_h12_inf')
        f = ina_h12_inf.pyfunc()

    # Calculate inactivation curve
    vs = np.linspace(-160, 40, 200)
    ss = f(vs)

    # Perform simulations
    v1s, v2s = [], []
    a1s, a2s = [], []
    for k in ko_levels:
        vr = state_fast(model, protocol, k)[vm]
        v1s.append(vr)
        a1s.append(f(vr))

        vr = state_slow(model, protocol, k)[vm]
        v2s.append(vr)
        a2s.append(f(vr))

    # Done!
    return vs, ss, v1s, a1s, v2s, a2s


def limit_cycle(model, protocol, rel_tol=1e-5, max_beats=4000, max_period=10):
    """
    Pre-paces a model to periodic orbit ("steady state").

    Arguments
    ``model``
        The model to pre-pace
    ``protocol``
        The protocol to pre-pace with
    ``rel_tol``
        The tolerance with which to judge a state as steady
    ``max_beats``
        The maximum number of beats to try
    ``max_period``
        The number of successive APs to check against. Alternans with
        `period < max_period` should be detected by this method.

    Return a tuple ``(state, period)`` where ``period`` is an integer
    indicating the periodicity of the states; 0 for unstable, 1 for stable, and
    2 or higher for alternans.

    Does not modify the model!
    """
    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    cl = protocol.characteristic_time()

    # Get scale of each state
    states = list(model.states())
    d = s.run(cl, log=myokit.LOG_STATE)
    x = np.array([d[var] for var in states])
    scale = np.max(x, axis=1) - np.min(x, axis=1)
    scale[scale == 0] = 1

    # Check if already at steady-state
    d = s.run(2 * cl, log_interval=cl, log=myokit.LOG_STATE)
    x = np.array([d[var] for var in states]).T
    dx = np.abs(x[0] - x[1]) / scale
    if np.max(dx) < rel_tol:
        return s.state(), 1

    beats = 0
    period = 0
    duration = max_period * cl
    while beats < max_beats:

        # Run and capture a number of beats
        d = s.run(duration, log_interval=cl, log=myokit.LOG_STATE)
        beats += max_period
        x = np.array([d[var] for var in states]).T

        # Check to see if any of the beats are a repeat of the first beat
        period = 0
        for i in range(1, max_period):
            dx = np.abs(x[0] - x[i]) / scale
            if np.max(dx) < rel_tol:
                period = i
                break

        # Repeat found! Check if there's a second repeat
        if period > 0:
            # Simulate more beats if needed
            if 2 * period >= max_period:
                d = s.run(duration, log_interval=cl, log=myokit.LOG_STATE)
                beats += max_period
                x = np.concatenate((x, np.array([d[var] for var in states]).T))

            dx = np.abs(x[period] - x[2 * period]) / scale
            if np.max(dx) < rel_tol:
                print('Terminating after ' + str(beats) + ' beats')
                break
            else:
                period = 0

    if period > 1:
        print('WARNING: Detected alternans with period ' + str(period) + '.')
    elif period == 0:
        print('WARNING: Terminating after maximum number of beats.')
        dx = np.abs(x[0] - x[i]) / scale
        print('Final dx: ' + str(np.max(dx)))

    return s.state(), period


def log_fast(model, protocol, ko, lw=False):
    """
    Returns a simulation log starting 1 beat after a change to the given Ko.

    Set ``lw`` to True to only gather a state, no logs.
    """
    path, cl = _path_cl(model, protocol, ko)
    path1 = path + '-fast-state.zip'
    path2 = path + '-fast-log.zip'
    if lw:
        if os.path.isfile(path1):
            return
    else:
        if os.path.isfile(path1) and os.path.isfile(path2):
            return myokit.DataLog.load(path2).npview()

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)

    # Change Ko
    s.set_constant(model.labelx('K_o'), ko)

    # Pre-pace for 1 beat
    s.run(cl, log=myokit.LOG_NONE)

    # Store state
    state = s.state()
    myokit.save_state_bin(path1, state)
    if lw:
        return

    # Run for 1 second and save
    #d = s.run(cl)
    d = s.run(800)
    d.save(path2)

    return d.npview()


def log_slow(model, protocol, ko, lw=False):
    """
    Returns a simulation log starting 900 beats after a change to the given Ko.

    Set ``lw`` to True to only gather a state, no logs.
    """
    path, cl = _path_cl(model, protocol, ko)
    path1 = path + '-slow-state.zip'
    path2 = path + '-slow-log.zip'
    path3 = path + '-prog-log.zip'
    if os.path.isfile(path1):
        if lw:
            return
        elif os.path.isfile(path2) and os.path.isfile(path3):
            return myokit.DataLog.load(path2).npview()

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)

    # Change Ko
    s.set_constant(model.labelx('K_o'), ko)

    # Pre-pace for 900 beats
    log = myokit.LOG_NONE if lw else None
    e = s.run(900 * cl, log=log, log_interval=cl)

    # Store progress log
    if not lw:
        e.save(path3)
    del(e)

    # Store state
    state = s.state()
    myokit.save_state_bin(path1, state)
    if lw:
        return

    # Run for 1 beat and save
    #d = s.run(cl)
    d = s.run(800)
    d.save(path2)

    return d.npview()


def log_progress(model, protocol, ko):
    """
    Returns a log showing the progress of a simulation over 900 beats, where
    only the first time point of each beat is logged.
    """
    path, cl = _path_cl(model, protocol, ko)
    path += '-prog-log.zip'

    if not os.path.isfile(path):
        log_slow(model, protocol, ko)

    # Return log
    return myokit.DataLog.load(path).npview()


def model(name=None):
    """
    Returns a Myokit model, depending on the 1st command line argument to the
    script.
    """
    if name is not None:
        return myokit.load_model(os.path.join('models', model_names[name]))

    try:
        model = model_names[sys.argv[1]]
    except KeyError:
        print('Unknown model: ' + sys.argv[1])
        print('Choose one of: ' + ', '.join([str(x) for x in model_names]))
        sys.exit(1)
    except IndexError:
        print(f'Syntax: {sys.argv[0]} model')
        print('Choose one of: ' + ', '.join([str(x) for x in model_names]))
        sys.exit(1)

    return myokit.load_model(os.path.join('models', model))


def models():
    """Returns a mapping from names to myokit.Model instances."""
    models = {}
    for name, fname in model_names.items():
        models[name] = myokit.load_model(os.path.join('models', fname))
    return models


def _path_cl(model, protocol, ko=None):
    """
    Returns a path and cycle length for the given model, protocol, and Ko.
    """

    # Get info from model
    name = model.name()

    # Get info from protocol
    if len(protocol) != 1:
        raise ValueError('Expecting a protocol with a single event')
    e = protocol.head()
    cl = e.period()
    offset = e.start()
    duration = e.duration()
    level = e.level()

    # Create filename
    path = f'{name}-cl-{fstr(cl)}-o-{fstr(offset)}-d-{fstr(duration)}'
    if level != 1:
        path += '-level-' + fstr(level)
    if ko is not None:
        path += '-ko-' + fstr(ko)
    path = os.path.join(sdata, path)

    return path, cl


def prepare_model(
        model, protocol=None, fix_cleft_ko=True, pre_pace=True, pA=False):
    """
    Prepares a model by setting the desired units, adding a voltage-clamp
    switch, and pre-pacing.

    The following variables / labels guaranteed to exist after this:

    - A time variable in ms
    - ``membrane_potential`` in mV
    - ``membrane_capacitance`` (units unchanged)
    - ``K_i`` in mM
    - ``K_o`` in mM
    - ``Na_i`` in mM
    - ``cellular_current`` in A/F
    - ``cellular_current_nA`` in nA (for R in mV/nA = MOhm)
    - ``stimulus_current`` in A/F
    - ``diffusion_current`` in A/F
    - ``I_CaL`` in A/F
    - ``I_K1`` in A/F
    - ``I_Kr`` in A/F
    - ``I_Kur`` in A/F
    - ``I_NaK`` in A/F
    - ``I_to`` in A/F

    And the following will be created
    - ``clamp`` (dimensionless)

    If ``fix_cleft_ko=True`` the method will detect any cleft potassium state
    variable and fix it to the bulk external value.

    Pre-pacing can be disabled by setting ``pre_pace=False``.

    Currents can be converted to ``pA`` instead by setting ``pA=True``.
    """
    # Test if already called
    if model.label('clamp') is not None:
        print('Prepare_model already called')
        return

    # Create protocol
    if protocol is None:
        protocol = myokit.pacing.blocktrain(1000, duration=0.5, offset=50)

    # Get model variables
    t = model.timex()
    v = model.labelx('membrane_potential')
    C = model.labelx('membrane_capacitance')
    ki = model.labelx('K_i')
    ko = model.labelx('K_o')
    nai = model.labelx('Na_i')
    i_ion = model.labelx('cellular_current')
    i_stim = model.labelx('stimulus_current')
    i_diff = model.bindingx('diffusion_current')
    clabels = ['I_CaL', 'I_K1', 'I_Kr', 'I_Kur', 'I_NaK', 'I_to']
    cvars = [model.labelx(c) for c in clabels]

    # Check variables have units
    variables = cvars + [v, t, C, ki, ko, nai, i_ion, i_stim, i_diff]
    for var in variables:
        if var.unit() is None:
            raise ValueError(
                'No unit set for ' + str(var) + ' in ' + str(model))

    # Convert remaining currents
    i_unit = 'pA' if pA else 'A/F'
    helpers = [C.rhs()]
    v.convert_unit('mV')
    ki.convert_unit('mM')
    ko.convert_unit('mM')
    nai.convert_unit('mM')
    i_ion.convert_unit(i_unit, helpers=helpers)
    i_stim.convert_unit(i_unit, helpers=helpers)
    i_diff.convert_unit(i_unit, helpers=helpers)
    for qname in currents(model):
        var = model.get(qname)
        if var.unit() is None:
            raise ValueError(
                'No unit set for ' + str(var) + ' in ' + str(model))
        var.convert_unit(i_unit, helpers=helpers)
    t.convert_unit('ms')

    # Get version of i_ion in nA (so we get mV/nA = MOhm)
    i_ion_nA = i_ion.parent().add_variable_allow_renaming('i_ion_nA')
    i_ion_nA.set_rhs(myokit.Name(i_ion))
    i_ion_nA.set_unit(i_ion.unit())
    i_ion_nA.set_label('cellular_current_nA')
    i_ion_nA.convert_unit('nA', helpers=helpers)

    # Set baseline Ko to 5.4
    ko.set_rhs('5.4 [mM]')

    # Update the model with a voltage clamping switch
    clamp = v.add_variable('clamp')
    clamp.set_rhs(0)
    clamp.set_label('clamp')
    v.set_rhs('if(clamp, 0, ' + str(v.rhs()) + ')')

    # Fix cleft potassium to bulk level
    cache_suffix = None
    kc = model.label('K_c')
    if fix_cleft_ko and kc is not None:
        print('Fixing cleft-K to bulk external level: ' + model.name())
        # Fix to bulk level
        kc.demote()
        kc.set_rhs(myokit.Name(ko))
        # Move ko and kc into their own component, to reduce chance of this
        # causing component cycles (for the opencl sim)
        c = model.add_component_allow_renaming('external_potassium')
        kc.parent().move_variable(kc, c)
        ko.parent().move_variable(ko, c)
        cache_suffix = '-fixed-cleft'

    # Pre-pace
    rel_tol = 0.6 if 'koiv' in model.name() else 1e-5
    if pre_pace and 'koiv' not in model.name():
        print('Pre-pacing: ' + model.name())
        print('  with ' + str(protocol))
        model.set_state(cached_limit_cycle(
            model, protocol, suffix=cache_suffix, rel_tol=rel_tol)[0])


def quasi_instantaneous_iv(
        model, protocol, times, current_label='cellular_current',
        simulation=None):
    """
    Plots a "quasi instantaneous IV curve" in the style of [1], where all
    states are fixed at their current values, but INa activation is set to its
    steady-state for any given V.

    Returns a tuple ``(vts, its, iv_voltages, iv_currents)`` where ``vts`` is a
    list of membrane potentials reached at the requested ``times`` and ``its``
    is a list of the corresponding currents. The array ``iv_voltages`` contains
    the membrane potentials at which IV curves were evaluated, and
    ``iv_currents`` is a list of arrays with the corresponding current values.

    [1] Fink, Giles, Noble (2006) Contributions of inwardly rectifying K+
        currents to repolarization assessed using mathematical models of human
        ventricular myocytes
    """
    # Get simulation
    if simulation is None:
        simulation = myokit.Simulation(model, protocol)
        simulation.set_tolerance(1e-8, 1e-8)

    # Get clone of model, get variables
    fixed_model = myokit.lib.hh.convert_hh_states_to_inf_tau_form(model)
    fixed_states = list(fixed_model.states())
    fixed_v = fixed_model.labelx('membrane_potential')
    fixed_i = fixed_model.labelx(current_label)
    fixed_m = fixed_model.labelx('ina_m')
    try:
        finf = fixed_m['inf']
        ftau = fixed_m['tau']
    except KeyError:
        finf, ftau = myokit.lib.hh.get_inf_and_tau(fixed_m)

    # Demote all variables except v
    for var in fixed_states:
        if var != fixed_v:
            demote(var)

    # Unbind all variables except time
    for label, var in list(fixed_model.bindings()):
        if label != 'time':
            var.set_binding(None)

    # Set ina.m to equal its inf variable
    fixed_m.set_rhs(myokit.Name(finf))

    # Create iv curves
    iv_voltages = np.linspace(-110, 50)
    iv_currents = []
    vts = []
    its = []

    for time in times:
        # Simulate up to t
        simulation.reset()
        simulation.run(time, log=myokit.LOG_NONE)

        # Update fixed model with new state
        state = simulation.state()
        for i, var in enumerate(fixed_states):
            if var == fixed_v:
                var.set_state_value(state[i])
                vts.append(state[i])
            else:
                var.set_rhs(state[i])

        # Calculate IV
        f = fixed_i.pyfunc()
        iv_currents.append(f(iv_voltages))
        its.append(f(fixed_v.state_value()))

    # Return
    return vts, its, iv_voltages, iv_currents


def rm(model, protocol, dt, dv, tmax=1000, simulation=None):
    """
    Calculate the series resistance using the given ``model`` and ``protocol``,
    with a step duration of ``dt`` and magnitude of ``dv``.

    At each time ``t``, the membrane is clamped to ``V(t) + dv`` for ``dt``
    seconds. The procedure is then repeated for ``V(t) - dv``.

    Returns a tuple ``(time, resistance)``.
    """
    if simulation is None:
        s = myokit.Simulation(model, protocol)
        s.set_tolerance(1e-8, 1e-8)
    else:
        s = simulation
        s.reset()

    v = model.labelx('membrane_potential')
    i_ion = model.labelx('cellular_current_nA')
    clamp = model.labelx('clamp')

    ts = np.linspace(0, tmax, 1 + tmax // 4)

    rs = []
    for t in ts:
        s.reset()
        s.set_constant(clamp, 0)
        s.run(t, log=myokit.LOG_NONE)
        x = s.state()
        s.set_constant(clamp, 1)

        # Run with +dv for dt time units
        x1 = list(x)
        x1[v.indice()] += dv
        s.set_state(x1)

        d = s.run(dt, log=[i_ion]).npview()
        i1 = d[i_ion]

        # Run with -dv for dt time units
        x2 = list(x)
        x2[v.indice()] -= dv
        s.set_state(x2)
        s.set_time(t)
        d = s.run(dt, log=[i_ion]).npview()
        i2 = d[i_ion]

        # store R
        r = 2 * dv / (i1[-1] - i2[-1])
        rs.append(r)
    s.set_constant(clamp, 0)

    return np.array(ts), np.array(rs)


def _sd_log(model, level, long=False):
    """ Run an SD experiment. """

    # Uses default protocol, for consistency with other figures

    # Number of beats to wait after Ko change
    n_beats = 900 if long else 1
    printing = False

    # Try returning cached file
    path = os.path.join(
        sdata, f'sd-{model.name()}-ko-{fstr(level)}-b-{n_beats}.csv')
    try:
        print(f'Loaded data from {path}')
        return myokit.DataLog.load_csv(path)
    except Exception:
        pass

    # Prepare model
    prepare_model(model)
    protocol = default_protocol()
    cl = 1000

    # Get typical stimulus amplitude
    temp = model.clone()
    x = temp.binding('pace')
    x.set_binding(None)
    x.set_rhs(1)
    default_amplitude = -temp.labelx('stimulus_current').eval()
    print(f'Default amplitude: -{default_amplitude} A/F')
    del(temp, x)

    # Get some variables
    tt = model.time()
    vm = model.labelx('membrane_potential')
    st = model.labelx('stimulus_current')
    ko = model.labelx('K_o')

    # Update i_stim formulation
    sa = st.parent().add_variable_allow_renaming('amp')
    sa.set_unit('A/F')
    sa.set_rhs(default_amplitude)
    st.set_rhs(myokit.Multiply(
        myokit.PrefixMinus(myokit.Name(sa)),
        myokit.Name(model.bindingx('pace'))
    ))

    # Durations to test
    #durations = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    durations = np.logspace(-2, np.log10(5), 100)

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-12, 1e-14)
    s.set_constant(ko, level)

    # Pre-pace for 1s
    print(f'Pacing for {n_beats} beats at ko={level}')
    s.pre(n_beats * cl)

    # Now pre-pace until just before the stimulus
    s.pre(49)

    # Get an expected Vmax and APD
    d = s.run(cl - 50, log=[tt, vm]).npview()
    vmin = np.min(d[vm])
    vmax = np.max(d[vm])
    vacc = -30  # Minimum accepted max V
    vap = 0.9 * vmin

    try:
        apd = d.apd(threshold=vap)['duration'][0]
    except IndexError:
        apd = float('nan')
    tmax = min(2 * apd, cl - 50)
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
        d = s.run(20, log=[tt, vm])
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

            amplitudes.append(amax)

    # Store
    d = myokit.DataLog()
    d['duration'] = durations
    d['amplitude'] = amplitudes
    d.save_csv(path)

    return d


def sd_log_fast(model, ko):
    """
    Returns a DataLog containing entries ``duration`` and ``amplitude`` that
    provide stimulus-duration data for the given model, one beat after changing
    to the given ko.
    """
    return _sd_log(model, ko, False)


def sd_log_slow(model, ko):
    """
    Returns a DataLog containing entries ``duration`` and ``amplitude`` that
    provide stimulus-duration data for the given model, 900 beats after
    changing to the given ko.
    """
    return _sd_log(model, ko, True)


def splash(figure_name, model_name=None):
    """Prints a figure name to the screen, in a flashy way."""
    if model_name is not None:
        figure_name = f'{figure_name} with {model_name}'
    print()
    print('=' * 79)
    print(f' Generating {figure_name}')
    print('=' * 79)


def state_fast(model, protocol, ko, lw=False):
    """
    Returns the model state 1 beat after a change to the given [K]o.

    Set ``lw`` to True to run a lightweight simulation, if needed, that doesn't
    gather or store simulation logs.
    """
    path, cl = _path_cl(model, protocol, ko)
    path += '-fast-state.zip'

    # Simulate
    if not os.path.isfile(path):
        log_fast(model, protocol, ko, lw)

    # Return state
    return model.map_to_state(myokit.load_state_bin(path))


def state_slow(model, protocol, ko, lw=False):
    """
    Returns the model state 900 beats after a change to the given [K]o.

    Set ``lw`` to True to run a lightweight simulation, if needed, that doesn't
    gather or store simulation logs.
    """
    path, cl = _path_cl(model, protocol, ko)
    path += '-slow-state.zip'

    # Simulate
    if not os.path.isfile(path):
        log_slow(model, protocol, ko, lw)

    # Return state
    return model.map_to_state(myokit.load_state_bin(path))


def strand_conductance_and_stim_amplitude(model):
    """
    Returns a tuple ``g, f`` where ``g`` is the cell-to-cell conductance to
    obtain a CV of 80cm/s, while ``f`` is the minimum stimulus amplitude (1s
    stimulus) for that conductance.
    """

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
        raise NotImplementedError(f'No conductance known for {model.name()}.')

    return g, f


def strand_test_beat(model, conductance=True):
    """
    Performs a test beat, and shows the resulting CV.

    Arguments:

    ``model``
        The model to use
    ``conductance``
        ``True`` (default) to try a conductance, ``False`` to try a stimulus
        threshold instead.
    """
    # Prepare model, pre-pace with default protocol (not the one used here!)
    prepare_model(model)

    # Get conductance and stimulus level
    g, f = strand_conductance_and_stim_amplitude(model)

    # Set stimulus amplitude as multiple of minimum
    if conductance:
        factor = 2
        print(f'Using stimulus that is {factor} times the minimum amplitude.')
        f *= factor
        del(factor)

    # Create protocol
    protocol = myokit.pacing.blocktrain(
        period=1000, offset=50, duration=5, level=f)

    # Pre-pace without caching, but now using real protocol
    rel_tol = 0.6 if 'koiv' in model.name() else 1e-5
    print('Pre-pacing single cell model with new protocol')
    model.set_state(limit_cycle(model, protocol, rel_tol=rel_tol)[0])

    # Get some variables
    tt = model.time()
    vm = model.labelx('membrane_potential')

    # Create simulation
    s = myokit.SimulationOpenCL(
        model, protocol, ncells=(200), rl=True,
        precision=myokit.DOUBLE_PRECISION)
    s.set_step_size(0.001)
    s.set_conductance(g)

    # Test stimulus size and/or conductance
    print('Simulating test beat for '
          + ('conductance' if conductance else 'stimulus amplitude'))
    print(f'f={f}, g={g}')
    d = s.run(100, log=[tt, vm])
    print('Done')
    b = d.block1d()
    cv = b.cv(name=vm.qname(), length=0.01)
    print(f'CV: {cv} cm/s')


def strand_test_stability(model):
    """
    Checks the stability of a strand simulation after 50 beats pre-pacing.
    """
    # Number of strand pre-pacing beats
    n_pre = 50

    # Load cached state
    dpath = 'strand-data'
    fname = f'{model.name()}-state-{n_pre}.zip'
    path = os.path.join(dpath, fname)
    if not os.path.isfile(path):
        raise ValueError(f'No pre-paced state found. Missing {path}')
    print(f'Loading state after {n_pre} beats at baseline Ko')
    state_pre = myokit.load_state_bin(path)

    fname = f'{model.name()}-state-{n_pre + 1}.zip'
    path = os.path.join(dpath, fname)
    if os.path.isfile(path):
        print(f'Loading state after {n_pre + 1} beats at baseline Ko')
        state_post = myokit.load_state_bin(path)
    else:
        print(f'Simulating a beat at baseline Ko...')

        # Create protocol
        g, f = strand_conductance_and_stim_amplitude(model)
        protocol = myokit.pacing.blocktrain(
            period=1000, level=2 * f, duration=5, offset=50)

        # Prepare model, without pre-pacing
        prepare_model(model, protocol, pre_pace=False)

        # Create simulation
        s = myokit.SimulationOpenCL(
            model, protocol, ncells=(200), rl=True,
            precision=myokit.DOUBLE_PRECISION)
        s.set_step_size(0.001)
        s.set_conductance(g)
        s.set_state(state_pre)

        # Run
        s.run(1000, log=myokit.LOG_NONE)
        state_post = s.state()
        myokit.save_state_bin(path, state_post)

    state_pre = np.array(state_pre)
    state_post = np.array(state_post)
    e = np.abs((state_post - state_pre) / state_pre) * 100
    e = e[np.isfinite(e)]
    print(f'{np.round(np.max(e), 3)} %')


def strand_ap_fast(model, ko):
    """
    Returns a log of the 2nd beat after a Ko change, in a strand simulation.
    """
    cvt = 40
    i_final = 52
    path = os.path.join(
        'strand-data',
        f'{model.name()}-cvt-{cvt}-ko-{fstr(ko)}-beat{i_final}-cv.bin'
    )

    if not os.path.isfile(path):
        _strand_sim(model, ko)
    return myokit.DataBlock1d.load(path).to_log().npview()


def strand_ap_slow(model, ko):
    """
    Returns a log of the 901st beat after a Ko change, in a strand simulation.
    """
    cvt = 40
    i_final = 951
    path = os.path.join(
        'strand-data',
        f'{model.name()}-cvt-{cvt}-ko-{fstr(ko)}-beat{i_final}-cv.bin')

    if not os.path.isfile(path):
        _strand_sim(model, ko)
    return myokit.DataBlock1d.load(path).to_log().npview()


def strand_stats_fast(model, ko):
    """
    Returns a log of derived stats for the 2nd beat after a Ko change, in a
    strand simulation.
    """
    cvt = 40
    i_final = 52
    path = os.path.join(
        'strand-data',
        f'{model.name()}-cvt-{cvt}-ko-{fstr(ko)}-beat{i_final}-wavelength.csv'
    )

    if not os.path.isfile(path):
        _strand_sim(model, ko)
    return myokit.DataLog.load_csv(path)


def strand_stats_slow(model, ko):
    """
    Returns a log of derived stats for the 951st beat after a Ko change, in a
    strand simulation.
    """
    cvt = 40
    i_final = 951
    path = os.path.join(
        'strand-data',
        f'{model.name()}-cvt-{cvt}-ko-{fstr(ko)}-beat{i_final}-wavelength.csv'
    )

    if not os.path.isfile(path):
        _strand_sim(model, ko)
    return myokit.DataLog.load_csv(path)


def _strand_sim(model, level, cvt=40, long=False):

    # Get path for data
    dpath = 'strand-data'
    name = f'{model.name()}-cvt-{cvt}-ko-{fstr(level)}'

    # Show details
    print('=' * 79)
    print(f'Running for {model.name()} with [K]o-change to {level}')
    print(f'Propagating wave threshold: CV > {cvt}')
    print(f'Running {"LOOOOONG 15min" if long else "1s"} simulation')

    # Stimulus protocol timing
    cl = 1000
    t0 = 50
    sd = 5
    print(f'Stimulus of {sd} ms applied every {cl} ms with a {t0} offset.')

    # Number of strand pre-pacing beats
    n_pre = 50

    # Conduction velocity settings
    cv_cell_length = 0.01   # In cm

    #
    # Note: The code below is almost independent of the choice of variables
    # above.
    # However: It assumes that n_pre is constant, and refers to simulations
    # with some number of n_pre and some number of n_beats as
    # sim_{n_pre+n_beats}. But if n_pre varies, this becomes blurry.
    # TODO: Either re-write as having n_pre fixed, or change the file format to
    # differentiate between n_pre and n_beats.
    #

    # Set conductance and stimulus level
    g, f = strand_conductance_and_stim_amplitude(model)

    # Create protocol (twice the minimum stimulus size)
    pargs = dict(period=cl, level=2 * f, duration=sd)
    protocol = myokit.pacing.blocktrain(offset=t0, **pargs)

    # Prepare model, without pre-pacing
    prepare_model(model, protocol, pre_pace=False)

    # Pre-pace, check if cached result was used
    rel_tol = 0.6 if 'koiv' in model.name() else 1e-5
    state, use_cached = cached_limit_cycle(model, protocol, rel_tol=rel_tol)

    # Get some variables
    tt = model.time()
    vm = model.labelx('membrane_potential')
    ko = model.labelx('K_o')

    # CV method arguments
    cvargs = dict(name=vm.qname(), length=cv_cell_length)

    # Create simulation
    s = myokit.SimulationOpenCL(
        model, protocol, ncells=(200), rl=True,
        precision=myokit.DOUBLE_PRECISION)
    s.set_step_size(0.001)
    s.set_conductance(g)

    # Pre-pace for n_pre beats at baseline ko
    fname = f'{model.name()}-state-{n_pre}.zip'
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
    fname1 = f'{name}-state{i_ready}.zip'
    fname2 = f'{name}-beat{i_ready}.bin'
    path1 = os.path.join(dpath, fname1)
    path2 = os.path.join(dpath, fname2)
    use_cached &= os.path.isfile(path1)
    use_cached &= os.path.isfile(path2)
    print(path1, os.path.isfile(path1), use_cached)
    print(path2, os.path.isfile(path2), use_cached)
    if use_cached:
        print(f'Loading state {i_ready} (state after {n_beats} beat(s) at new'
              f' Ko)')
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
    if long:
        n_done = n_beats
        n_beats = 15 * 60           # 15 minutes, assuming CL=1000ms
        i_ready = n_pre + n_beats   # Index of pre-paced beat
        i_final = i_ready + 1       # Index of test beat
        n_todo = n_beats - n_done

        fname1 = f'{name}-state{i_ready}.zip'
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
                i_hi = (n_pre + n_todo) // 50 - 1
                for i in range(i_hi, i_lo, -1):
                    j = 50 * i
                    fname2 = f'{name}-state{j}.zip'
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
            # At this point, we have loaded an intermediate state, but still
            # need to simulate some more, so use_cached for future files should
            # be false.
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
                    fname2 = f'{name}-state{j}.zip'
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
    fname = f'{name}-beat{i_final}-cv.bin'
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
    fname1 = f'{name}-beat{i_final}-refractory-period.csv'
    fname2 = f'{name}-beat{i_final}-wavelength.csv'
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
            s.set_protocol(
                myokit.pacing.blocktrain(offset=t, limit=1, **pargs))

            # Run until a short while after t
            # At 5cm/s, it would take 1.5cm/(5cm/s) = 0.3s for cell 150 to get
            # activated, so 300ms sim should be long enough for even very slow
            # CVs to be measured (cv method ignores first and last 50 cells).
            # At 10cm/s, 150ms should be enough
            t_final = t + 150
            duration = t_final - t0
            print(f'  Simulating from t={t0} ms to t={t_final} ms for a beat'
                  f' at t={t} ms (duration {duration} ms).')
            d = s.run(duration, log=[tt, vm])
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
                # RP = Lowest conducting t, minus t start, minus stimulus
                # offset
                rp = t2 - time_ready - t0
                wl = cv * rp * 0.01  # In mm

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
        print(f'No FRP known: NOT simulating beat {i_final} and refractory'
              f' beat.')
    else:
        fname = f'{name}-beat{i_final}-refrac.bin'
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
            s.set_protocol(
                myokit.pacing.blocktrain(offset=t, limit=1, **pargs))
            duration = t_final - tmin
            print(f'Simulating from t={tmin} ms to t={t_final} ms for a beat'
                  f' at t={t} ms (duration {duration} ms).')
            d = s.run(duration, log=d)

            print('Storing final and refractory beat')
            b_final_refrac = d.block1d()
            b_final_refrac.save(path)
            del(d)
        del(fname, path)


def strand_frp_graph(ax, model, slow=False):
    """
    Create a plot of the refractory period on the given axes.
    """
    vm = model.labelx('membrane_potential')

    cl = 1000
    offset = 50

    # Update axes
    ax.minorticks_on()
    if 'voigt' in model.name():
        ax.set_ylim(-85, 35)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('V cell 1 (mV)')

    for i, level in enumerate(ko_levels):
        print(f'Creating plot for ko={level}')

        # Load data
        if slow:
            d = strand_ap_slow(model, level)
            rp = strand_stats_slow(model, level)['rp'][0]
        else:
            d = strand_ap_fast(model, level)
            rp = strand_stats_fast(model, level)['rp'][0]

        # Plot AP
        ax.plot(
            d.time() * 1e-3, d[vm, 0], color=ko_colors[i], label=f'{level} mM')

        # Plot FRP
        if np.isfinite(rp):
            frp = (d.time()[0] + offset + rp) * 1e-3
            ax.axvline(frp, color=ko_colors[i], alpha=0.5)

    # Finalise plot
    ticks = [int(d.time()[0] * 1e-3)] * 3
    ticks[1] += cl * 1e-3 / 2
    ticks[2] += int(cl * 1e-3)

    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) for x in ticks])
    ax.legend(loc='upper right')


def strand_stats_graphs(
        ax1, ax2, ax3, ax4, ax5, model, long=False, crosses=[]):
    """
    Plots overview graphs of apd, rp, cv, wl, and Vr on the given axes.

    Concentrations at which to show crosses instead of squares can be passed in
    with ``crosses``

    """
    # Gather data
    vm = model.labelx('membrane_potential')
    ap, rp, cv, wl, vr = [], [], [], [], []
    for i, level in enumerate(ko_levels):

        # Load data
        if long:
            d = strand_ap_slow(model, level)
            csv = strand_stats_slow(model, level)
        else:
            d = strand_ap_fast(model, level)
            csv = strand_stats_fast(model, level)

        # Gather stats
        ap.append(csv['apd'][0])  # ms
        rp.append(csv['rp'][0])   # ms
        cv.append(csv['cv'][0])   # cm/s
        wl.append(csv['wl'][0])   # mm
        vr.append(d[vm, 0][0])

    # Shared limits
    xlim = 2.3, 5.6

    # Vr
    ax = ax1
    ax.set_ylabel('$V_r$ (mV)')
    ax.set_xlim(*xlim)
    if 'voigt' in model.name():
        ax.set_ylim(-78, -49)
    ax.set_xticks(ko_levels)
    ax.set_xticklabels([])
    ax.plot(ko_levels, vr, '-')
    crosses = set(crosses)
    for x, y, c in zip(ko_levels, vr, ko_colors):
        ax.plot(x, y, 'x' if x in crosses else 's', color=c)

    # CV
    ax = ax2
    ax.set_ylabel('CV (cm/s)')
    ax.set_xlim(*xlim)
    if 'voigt' in model.name():
        ax.set_ylim(-10, 86)
        ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_xticks(ko_levels)
    ax.set_xticklabels([])
    ax.plot(ko_levels, cv, '-')
    for x, y, c in zip(ko_levels, cv, ko_colors):
        ax.plot(x, y, 'x' if x in crosses else 's', color=c)

    # APD
    ax = ax3
    ax.set_ylabel('APD (ms)')
    ax.set_xlim(*xlim)
    if 'voigt' in model.name():
        ax.set_ylim(200, 600)
        ax.set_yticks([200, 300, 400, 500, 600])
    ax.set_xticks(ko_levels)
    ax.set_xticklabels([])
    ax.plot(ko_levels, rp, '-')
    for x, y, c in zip(ko_levels, rp, ko_colors):
        ax.plot(x, y, 'x' if x in crosses else 's', color=c)

    # FRP
    ax = ax4
    ax.set_ylabel('FRP (ms)')
    ax.set_xlim(*xlim)
    if 'voigt' in model.name():
        ax.set_ylim(200, 600)
        ax.set_yticks([200, 300, 400, 500, 600])
    ax.set_xticks(ko_levels)
    ax.set_xticklabels([])
    ax.plot(ko_levels, rp, '-')
    for x, y, c in zip(ko_levels, rp, ko_colors):
        ax.plot(x, y, 'x' if x in crosses else 's', color=c)

    # WL
    ax = ax5
    ax.set_xlabel('External potassium concentration $[K^+]_o$')
    ax.set_ylabel('WL (mm)')
    ax.set_xlim(*xlim)
    if 'voigt' in model.name():
        ax.set_ylim(160, 400)
        ax.set_yticks([200, 300, 400])
    ax.set_xticks(ko_levels)
    ax.plot(ko_levels, wl, '-')
    for x, y, c in zip(ko_levels, wl, ko_colors):
        ax.plot(x, y, 'x' if x in crosses else 's', color=c)

