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


# Models
model_names = {
    'courtemanche': 'courtemanche-1998.mmt',
    'grandi': 'grandi-2011.mmt',
    'koivumaki': 'koivumaki-2011.mmt',
    'maleckar': 'maleckar-2008.mmt',
    'ni': 'ni-2017.mmt',
    'nygren': 'nygren-1998.mmt',
    'voigt': 'voigt-heijman-2013.mmt',
    #'modified': 'modified-grandi-2011.mmt',
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


def figure_dir():
    """
    Returns the path to a directory for storing figures, depending on the
    script name.
    """
    path = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    if path[:1] == 'f':
        path = 'figure-' + path[1:]
    else:
        path = 'sup-figure-' + path
    if not os.path.isdir(path):
        os.makedirs(path)
    return path


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
        print(f'Choose one of: ' + ', '.join([str(x) for x in model_names]))
        sys.exit(1)
    except IndexError:
        print(f'Syntax: {sys.argv[0]} model')
        print(f'Choose one of: ' + ', '.join([str(x) for x in model_names]))
        sys.exit(1)

    return myokit.load_model(os.path.join('models', model))


def models():
    """Returns a mapping from names to myokit.Model instances."""
    models = {}
    for name, fname in model_names.items():
        models[name] = myokit.load_model(os.path.join('models', fname))
    return models


def prepare_model(
        model, protocol=None, fix_cleft_ko=False, pre_pace=True, pA=False):
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

    # Pre-pace
    if pre_pace and not 'koiv' in model.name():
        print('Pre-pacing: ' + model.name())
        model.set_state(limit_cycle(model, protocol))
        print(model.format_state(model.state()))
    else:
        print('NOT Pre-pacing: ' + model.name())


def limit_cycle(model, protocol, cl=None, rel_tol=1e-5, max_beats=4000,
                max_period=10, path=None):
    """
    Pre-paces a model to periodic orbit ("steady state").

    Arguments
    ``model``
    ``protocol``
    ``cl``
    ``rel_tol``
    ``max_beats``
    ``max_period``
    ``path``
    """

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)
    if cl is None:
        cl = protocol.characteristic_time()

    # Load steady-state from file, if given
    try:
        loaded = myokit.load_state(path)
        s.set_state(loaded)
    except Exception:
        loaded = None
        pass
    else:
        print('Loaded state from ' + str(path))

    # Get scale of each state
    states = list(model.states())
    d = s.run(cl, log=myokit.LOG_STATE)
    x = np.array([d[var] for var in states])
    scale = np.max(x, axis=1) - np.min(x, axis=1)
    scale[scale==0] = 1

    # Check if already at steady-state
    d = s.run(2 * cl, log_interval=cl, log=myokit.LOG_STATE)
    x = np.array([d[var] for var in states]).T
    dx = np.abs(x[0] - x[1]) / scale
    if np.max(dx) < rel_tol:
        return s.state() if loaded is None else loaded

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
            if 2 * period > max_period:
                d = s.run(duration, log_interval=cl, log=myokit.LOG_STATE)
                beats += max_period
                x = np.concatenate((x, np.array([d[var] for var in states]).T))

            dx = np.abs(x[period] - x[2 * period]) / scale
            if np.max(dx) < rel_tol:
                print('Terminating after ' + str(beats) + ' beats')
                break
            else:
                period = 0

    # Save state to file
    if path is not None:
        print('Saving final state to ' + str(path))
        myokit.save_state(path, s.state())

    if period > 1:
        print('WARNING: Detected alternans with period ' + str(period) + '.')
    elif period == 0:
        print('WARNING: Terminating after maximum number of beats.')
        dx = np.abs(x[0] - x[i]) / scale
        print('Final dx: ' + str(np.max(dx)))

    return s.state()


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

