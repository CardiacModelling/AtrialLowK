#
# Massive graphs of ionic fluxes and major affecting currents.
#
# Very model-specific.
#
import os
import matplotlib
import matplotlib.gridspec
import matplotlib.patches
import matplotlib.pyplot as plt
import myokit
import myokit.lib.plots
import numpy as np
import shared


def grandi(fpath, fname, model, pre_beats=0):
    """
    Ionic balance graph for Grandi model variants.
    """

    # Create protocol
    cl = 1000
    protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

    # Load model
    model = os.path.join('models', model)
    model = myokit.load_model(model)
    shared.prepare_model(model, protocol, pre_pace=False, pA=True)

    # Convert _jn  and _sl and such to pA too
    A_per_F = myokit.units.A / myokit.units.F
    C = model.labelx('membrane_capacitance')
    helpers = [C.rhs()]
    for var in model.variables(deep=True, state=False):
        if var.unit() == A_per_F and var.name()[-3:] in ('_jn', '_sl'):
            var.convert_unit('pA', helpers=helpers)

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)

    # Pre-pace?
    if pre_beats:
        print(f'Pre-pacing for {pre_beats} beats')
        s.pre(pre_beats * cl)

    # Create simulation: todo - use new simulation class
    print('Running simulation')
    d = s.run(1000).npview()

    # All currents, in pA
    currents = [
        'ikur.IKur',
        'ito.Ito',
        'ikr.IKr',
        'ical.ICaL_K',
        'ikp.IKp',
        'ik1.IK1',
        'stimulus.I_stim',
        'iks.IKs',

        'inak.INaK',

        'ina.INa',
        'inal.INaL',
        'ical.ICaL_Na',
        'inab.INaB',
        'inaca.INaCa',

        'ical.ICaL_Ca',
        'icab.ICaB',
        'ipca.IpCa',

        'iclca.IClCa',
        'iclb.IClB',
    ]
    ik = {
        'IKur': 1,
        'Ito': 1,
        'IKr': 1,
        'ICaL_K': 1,
        'IKp': 1,
        'IK1': 1,       # Positive = 1 K+ out
        'I_stim': 1,
        'INaK': -2,     # Positive = 2 K+ in, so should be times -2
        'IKs': 1,
    }
    ina = {
        'INa': 1,
        'INaK': 3,      # Positive = 3 Na+ out, so should be times 3
        'INaL': 1,
        'ICaL_Na': 1,
        'INaB': 1,      # Positive = 1 Na+ out
        'INaCa': 3,     # Positive = 3 Na+ out, so should be times 3
    }
    ica = {
        'ICaL_Ca': 0.5, # 1/2 an ion per unit of charge
        'ICaB': 0.5,
        'IpCa': 0.5,
        'INaCa': -1,    # Positive = 1 Ca2+ in, so should be times -1
    }
    icl = {
        'IClCa': -1,
        'IClB': -1,
    }

    labels = []
    dk = myokit.DataLog(time='time')
    dna = myokit.DataLog(time='time')
    dca = myokit.DataLog(time='time')
    dcl = myokit.DataLog(time='time')
    dk['time'] = d.time()
    dna['time'] = d.time()
    dca['time'] = d.time()
    dcl['time'] = d.time()

    # Calculate mol/s flux:
    #   F = 96485 C / mol
    #   A = C / s
    #   A/(zF) = mol / s
    #   pA/F = pmol / s
    #   pA/F = fmol / ms
    #   pA*1000/F = amol / ms
    #
    f = -1000 / 96485
    for current in currents:
        label = current.split('.')[1]
        labels.append(label)
        dk[label] = d[current] * ik.get(label, 0) * f
        dna[label] = d[current] * ina.get(label, 0) * f
        dca[label] = d[current] * ica.get(label, 0) * f
        dcl[label] = d[current] * icl.get(label, 0) * f

    # Colors for plotting
    cmap = plt.cm.get_cmap('tab20')
    colors = [cmap(i) for i in range(len(labels))]
    cmap = dict(zip(labels, colors))

    # Volume of bulk myoplasm in L
    Vmyo = model.get('geom.Vmyo').eval()
    Vsl = model.get('geom.Vsl').eval()
    Vjn = model.get('geom.Vjn').eval()
    Vsr = model.get('geom.Vsr').eval()

    # Total internal species
    ki = (
        + d['potassium.K_i'] * (Vmyo + Vsl + Vjn)       # mM/L * L = mmol
    ) * 1e15    # mmol * 1e15 = 1e-18 mol = attomol

    nai = (
        + d['sodium.Na_i'] * Vmyo
        + d['sodium.Na_sl'] * Vsl
        + d['sodium.Na_jn'] * Vjn
        + d['buffna.NaB_jn'] * Vjn
        + d['buffna.NaB_sl'] * Vsl
    ) * 1e15

    cai = (
        + d['calcium.Ca_i'] * Vmyo
        + d['calcium.Ca_sl'] * Vsl
        + d['calcium.Ca_jn'] * Vjn
        + d['calcium.Ca_sr'] * Vsr
        + d['calcium.Csqn'] * Vsr
        + d['buffca.TnCL'] * Vmyo
        + d['buffca.TnCHc'] * Vmyo
        + d['buffca.TnCHm'] * Vmyo
        + d['buffca.CaM'] * Vmyo
        + d['buffca.Myoc'] * Vmyo
        + d['buffca.Myom'] * Vmyo
        + d['buffca.SRB'] * Vmyo
        + d['buffca.SLL_sl'] * Vsl
        + d['buffca.SLH_sl'] * Vsl
        + d['buffca.SLL_jn'] * Vjn
        + d['buffca.SLH_jn'] * Vjn
    ) * 1e15

    cli = (
        + d['chloride.Cl_i'] * (Vmyo + Vsl + Vjn)
    ) * 1e15


    # Compare transmembrane flux with internal amount
    time = d.time()

    # K
    flux = (dk['INaK'][-1] + dk['IK1'][-1])
    deriv = (ki[-1] - ki[-2]) / (time[-1] - time[-2])
    print(f'K flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of K at t=1000: {deriv} amol/ms')

    # Na
    flux = (dna['INaK'][-1] + dna['INaCa'][-1] + dna['INaB'][-1])
    deriv = (nai[-1] - nai[-2]) / (time[-1] - time[-2])
    print(f'Na flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of Na at t=1000: {deriv} amol/ms')

    # Ca
    flux = (dca['IpCa'][-1] + dca['INaCa'][-1] + dca['ICaB'][-1])
    deriv = (cai[-1] - cai[-2]) / (time[-1] - time[-2])
    print(f'Ca flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of Ca at t=1000: {deriv} amol/ms')

    # Cl
    flux = +(dcl['IClB'][-1])
    deriv = (cli[-1] - cli[-2]) / (time[-1] - time[-2])
    print(f'Cl flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of Cl at t=1000: {deriv} amol/ms')



    #
    # Create figure
    #
    print('Plotting...')
    fig = plt.figure(figsize=(18, 11))
    fig.subplots_adjust(0.048, 0.045, 0.994, 0.945)  # Left, bottom, right, top
    grid = matplotlib.gridspec.GridSpec(12, 6, wspace=0.5, hspace=0.95)

    #
    # Relative contributions per species
    #
    ax00 = fig.add_subplot(grid[0:3, 0])
    myokit.lib.plots.cumulative_current(
        dk, labels, ax00, colors=colors, integrate=False, normalise=True)
    ax00.legend(loc=(0, 1.04), ncol=10).get_frame().set_alpha(1)
    ax00.set_ylabel('Contribution to inward flux')

    ax01 = fig.add_subplot(grid[0:3, 1])
    myokit.lib.plots.cumulative_current(
        dna, labels, ax01, colors=colors, integrate=False, normalise=True)

    ax02 = fig.add_subplot(grid[0:3, 2])
    myokit.lib.plots.cumulative_current(
        dca, labels, ax02, colors=colors, integrate=False, normalise=True)

    ax03 = fig.add_subplot(grid[0:3, 3])
    myokit.lib.plots.cumulative_current(
        dcl, labels, ax03, colors=colors, integrate=False, normalise=True)

    #
    # Absolute contributions: ions / s
    #
    ax10 = fig.add_subplot(grid[3:6, 0])
    ax10.set_ylabel('Selected currents\nInward K+ flux (amol/ms)')
    sel = ['IK1', 'INaK']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dk, sel, ax10, colors=col)

    ax11 = fig.add_subplot(grid[3:6, 1])
    ax11.set_ylabel('Selected currents\nInward Na+ flux (amol/ms)')
    sel = ['INaK', 'INaB', 'INaCa']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dna, sel, ax11, colors=col)

    ax12 = fig.add_subplot(grid[3:6, 2])
    ax12.set_ylabel('Selected currents\nInward Ca2+ flux (amol/ms)')
    sel = ['IpCa', 'INaCa', 'ICaB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dca, sel, ax12, colors=col)

    ax13 = fig.add_subplot(grid[3:6, 3])
    ax13.set_ylabel('Selected currents\nInward Cl- flux (amol/ms)')
    sel = ['IClB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dcl, sel, ax13, colors=col)

    #
    # Bulk amounts
    #
    def changeplot(ax, t, y):
        ax.axhline(y[0], color='#cccccc')
        ax.plot(t, y, color='k')

    ax20 = fig.add_subplot(grid[6:8, 0])
    ax20.set_ylabel('K+ ions (amol)')
    #ax20.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
    changeplot(ax20, d.time(), ki)

    ax21 = fig.add_subplot(grid[6:8, 1])
    ax21.set_ylabel('Na+ ions (amol)')
    changeplot(ax21, d.time(), nai)

    ax22 = fig.add_subplot(grid[6:8, 2])
    ax22.set_ylabel('Ca2+ ions (amol)')
    changeplot(ax22, d.time(), cai)

    ax23 = fig.add_subplot(grid[6:8, 3])
    ax23.set_ylabel('Cl- ions (amol)')
    changeplot(ax23, d.time(), cli)

    #
    # Background current
    #
    ax31 = fig.add_subplot(grid[8:10, 1])
    ax31.set_ylabel('INaB flux (amol/s)')
    sel = ['INaB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dna, sel, ax31, colors=col)

    ax32 = fig.add_subplot(grid[8:10, 2])
    ax32.set_ylabel('ICaB flux (amol/ms)')
    sel = ['ICaB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dca, sel, ax32, colors=col)

    ax33 = fig.add_subplot(grid[8:10, 3])
    ax33.set_ylabel('[Cl]i (mM)')
    ax33.set_ylabel('IClB flux (amol/ms)')
    sel = ['IClB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dcl, sel, ax33, colors=col)

    #
    # Phase-plots for concentrations
    #
    # Set `model` to final simulation state
    model.set_state(s.state())
    model.remove_derivative_references()

    def inf(y, x):
        m = model.clone()
        x = m.get(x)
        assert(x.is_state())
        m.binding('pace').set_binding(None)
        for state in list(m.states()):
            if state != x:
                state.set_rhs(state.state_value())
                state.demote()
        return m.get(y).pyfunc()

    color1 = 'blue'
    color2 = 'brown'

    f = inf('potassium.K_i', 'potassium.K_i')
    c = model.get('potassium.K_i').state_value()
    dc = 20
    x = np.linspace(c - dc, c + dc, 100)
    ax40 = fig.add_subplot(grid[10:12, 0])
    ax40.set_xlabel('[K+]i')
    ax40.set_ylabel('dot([K+]i)')
    ax40.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
    ax40.axvline(c, color='#cccccc')
    ax40.plot(x, f(x), color=color1)

    f1 = inf('sodium.Na_jn', 'sodium.Na_jn')
    c1 = model.get('sodium.Na_jn').state_value()
    f2 = inf('sodium.Na_sl', 'sodium.Na_sl')
    c2 = model.get('sodium.Na_sl').state_value()
    c = (c1 + c2) / 2
    dc = 0.001
    x = np.linspace(c - dc, c + dc, 100)
    ax41 = fig.add_subplot(grid[10:12, 1])
    ax41.set_xlabel('[Na+]')
    ax41.set_ylabel('dot([Na+])')
    ax41.axvline(c1, color='#cccccc')
    ax41.axvline(c2, color='#cccccc', ls='--')
    ax41.axhline(0, color='#777777')
    ax41.plot(x, f1(x), color=color1, label='jn')
    ax41.plot(x, f2(x), color=color2, label='sl', ls='--')
    ax41.legend()

    f1 = inf('calcium.Ca_jn', 'calcium.Ca_jn')
    c1 = model.get('calcium.Ca_jn').state_value()
    f2 = inf('calcium.Ca_sl', 'calcium.Ca_sl')
    c2 = model.get('calcium.Ca_sl').state_value()
    c = (c1 + c2) / 2
    dc = 0.6e-4
    x = np.linspace(max(0, c - dc), c + dc, 100)
    ax42 = fig.add_subplot(grid[10:12, 2])
    ax42.set_xlabel('[Ca2+]')
    ax42.set_ylabel('dot([Ca2+])')
    ax42.axvline(c1, color='#cccccc')
    ax42.axvline(c2, color='#cccccc', ls='--')
    ax42.axhline(0, color='#777777')
    ax42.plot(x, f1(x), color=color1, label='jn')
    ax42.plot(x, f2(x), color=color2, label='sl', ls='--')

    f = inf('chloride.Cl_i', 'chloride.Cl_i')
    c = model.get('chloride.Cl_i').state_value()
    dc = 5
    x = np.linspace(c - dc, c + dc, 100)
    ax43 = fig.add_subplot(grid[10:12, 3])
    ax43.set_xlabel('[Cl-]i')
    ax43.set_ylabel('dot([Cl-]i)')
    ax43.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
    ax43.axvline(c, color='#cccccc')
    ax43.plot(x, f(x), color=color1)

    #
    # Influence of concentrations on currents
    #
    sub = matplotlib.gridspec.GridSpecFromSubplotSpec(
        6, 3, subplot_spec=grid[:, 4:], wspace=0.5, hspace=0.5)

    def plot_inf(ax, y, x):

        dual = False
        if x == 'membrane.V':
            ax.set_xlabel('Vr (mV)')
            dx = 4
        elif x == 'potassium.K_i':
            ax.set_xlabel('[K+]i')
            dx = 10
        elif x == 'chloride.Cl_i':
            ax.set_xlabel('[Cl-]i')
            dx = 2
        elif x == 'sodium.Na':
            dual = True
            ax.set_xlabel('[Na+] jn,sl')
            dx = 1
        elif x == 'calcium.Ca':
            dual = True
            ax.set_xlabel('[Ca2+] jn,sl')
            dx = 0.0005
        else:
            raise NotImplementedError(x)

        label = y.split('.')[1]
        ax.set_ylabel(f'{label} (pA)')
        color = cmap[y.split('.')[1]]
        #f = -1000 / 96485
        #dk[label] = d[current] * ik.get(label, 0) * f

        if dual:
            x1, y1 = x + '_jn', y + '_jn'
            x2, y2 = x + '_sl', y + '_sl'
            c1 = model.get(x1).state_value()
            c2 = model.get(x2).state_value()
            xr = max(0, min(c1 - dx, c2 - dx)), max(c1 + dx, c2 + dx)
            r = np.linspace(xr[0], xr[1], 100)
            f1 = inf(y1, x1)
            f2 = inf(y2, x2)
            ax.axvline(c1, color='#cccccc')
            ax.axvline(c2, color='#cccccc', ls='--')
            ax.axhline(f1(c1), color='#cccccc')
            ax.axhline(f2(c2), color='#cccccc', ls='--')
            ax.plot(r, f1(r), color=color, label='jn', lw=2)
            ax.plot(r, f2(r), color=color, label='sl', lw=2, ls='--')
            yc = 0.5 * (f1(c1) + f2(c2))
        else:
            c = model.get(x).state_value()
            xr = c - dx, c + dx
            r = np.linspace(xr[0], xr[1], 100)
            f = inf(y, x)
            ax.axvline(c, color='#cccccc')
            ax.axhline(f(c), color='#cccccc')
            ax.plot(r, f(r), color=color, lw=3)
            yc = f(c)

        # Set scaling
        ax.set_xlim(*xr)
        dy = 10
        ax.set_ylim(yc - dy, yc + dy)


    # IK1(K, V)
    ax = fig.add_subplot(sub[0, 0])
    plot_inf(ax, 'ik1.IK1', 'potassium.K_i')
    ax.text(1.1, 0.5, 'Positive:\n1 K+ out', transform=ax.transAxes)
    ax = fig.add_subplot(sub[0, 2])
    plot_inf(ax, 'ik1.IK1', 'membrane.V')

    # INaK(Na, V)
    ax = fig.add_subplot(sub[1, 0])
    plot_inf(ax, 'inak.INaK', 'sodium.Na')
    ax.legend(loc=(1.1, 0))
    ax.text(1.1, 0.5, 'Positive:\n3 Na+ out, 2 K+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[1, 2])
    plot_inf(ax, 'inak.INaK', 'membrane.V')

    # INaB(Na, V)
    ax = fig.add_subplot(sub[2, 0])
    plot_inf(ax, 'inab.INaB', 'sodium.Na')
    ax.legend(loc=(1.1, 0))
    ax.text(1.1, 0.5, 'Negative: 1 Na+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[2, 2])
    plot_inf(ax, 'inab.INaB', 'membrane.V')

    # INaCa(Na, Ca, V)
    ax = fig.add_subplot(sub[3, 0])
    plot_inf(ax, 'inaca.INaCa', 'sodium.Na')
    ax = fig.add_subplot(sub[3, 1])
    plot_inf(ax, 'inaca.INaCa', 'calcium.Ca')
    ax.text(-0.45, 1.07, 'Negative: 1 Ca2+ out, 3 Na+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[3, 2])
    plot_inf(ax, 'inaca.INaCa', 'membrane.V')

    # ICaB(Ca, V)
    ax = fig.add_subplot(sub[4, 0])
    plot_inf(ax, 'icab.ICaB', 'calcium.Ca')
    ax.legend(loc=(1.1, 0))
    ax.text(1.1, 0.5, 'Negative: 1 Ca2+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[4, 2])
    plot_inf(ax, 'icab.ICaB', 'membrane.V')

    # IClB(Cl, V)
    ax = fig.add_subplot(sub[5, 0])
    plot_inf(ax, 'iclb.IClB', 'chloride.Cl_i')
    ax.text(1.1, 0.5, 'Negative: 1 Cl- out', transform=ax.transAxes)
    ax = fig.add_subplot(sub[5, 2])
    plot_inf(ax, 'iclb.IClB', 'membrane.V')

    #
    # Shading
    #
    # Add model name
    x = 0.672
    y = 0.182
    fig.patches.append(matplotlib.patches.Rectangle(
        (0, 0), 1, y, color='#eeeeee', fill=True, zorder=-1,
        transform=fig.transFigure))
    fig.patches.append(matplotlib.patches.Rectangle(
        (x, 0), 1 - x, 1, color='#eeeeee', fill=True, zorder=-1,
        transform=fig.transFigure))

    #
    # Store
    #
    path = os.path.join(fpath, fname)
    print(f'Writing to {path}')
    plt.savefig(path + '.png')


def ni(fpath, fname, model, pre_beats=0):
    """
    Ionic balance graph for Ni model variants.
    """

    # Create protocol
    cl = 1000
    protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

    # Load model
    model = os.path.join('models', model)
    model = myokit.load_model(model)
    shared.prepare_model(model, protocol, pre_pace=False, pA=True)

    # Convert _jn  and _sl and such to pA too
    A_per_F = myokit.units.A / myokit.units.F
    C = model.labelx('membrane_capacitance')
    helpers = [C.rhs()]
    for var in model.variables(deep=True, state=False):
        if var.unit() == A_per_F and var.name()[-3:] in ('_jn', '_sl'):
            var.convert_unit('pA', helpers=helpers)

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-8, 1e-8)

    # Pre-pace?
    if pre_beats:
        print(f'Pre-pacing for {pre_beats} beats')
        s.pre(pre_beats * cl)

    # Create simulation: todo - use new simulation class
    print('Running simulation')
    d = s.run(1000).npview()

    # All currents, in pA
    currents = [
        'ikur.IKur',
        'ito.Ito',
        'ikr.IKr',
        #'ikp.IKp',
        'ik1.IK1',
        'stimulus.i_stim',
        'iks.IKs',
        'inak.INaK',
        'ina.INa',
        'ibna.IbNa',
        'inaca.INaCa',
        'ical.ICaL',
        'ibca.IbCa',
        'icap.ICap',
    ]
    ik = {
        'IKur': 1,
        'Ito': 1,
        'IKr': 1,
        'IK1': 1,       # Positive = 1 K+ out
        'I_stim': 1,
        'IKs': 1,
        'INaK': -2,     # Positive = 2 K+ in, so should be times -2
    }
    ina = {
        'INaK': 3,      # Positive = 3 Na+ out, so should be times 3
        'INa': 1,
        'IbNa': 1,      # Positive = 1 Na+ out
        'INaCa': 3,     # Positive = 3 Na+ out, so should be times 3
    }
    ica = {
        'INaCa': -1,    # Positive = 1 Ca2+ in, so should be times -1
        'ICaL': 0.5, # 1/2 an ion per unit of charge
        'IbCa': 0.5,
        'ICap': 0.5,
    }

    labels = []
    dk = myokit.DataLog(time='time')
    dna = myokit.DataLog(time='time')
    dca = myokit.DataLog(time='time')
    dk['time'] = d.time()
    dna['time'] = d.time()
    dca['time'] = d.time()

    # Calculate mol/s flux:
    #   F = 96485 C / mol
    #   A = C / s
    #   A/(zF) = mol / s
    #   pA/F = pmol / s
    #   pA/F = fmol / ms
    #   pA*1000/F = amol / ms
    #
    f = -1000 / 96485
    for current in currents:
        label = current.split('.')[1]
        labels.append(label)
        dk[label] = d[current] * ik.get(label, 0) * f
        dna[label] = d[current] * ina.get(label, 0) * f
        dca[label] = d[current] * ica.get(label, 0) * f

    # Colors for plotting
    cmap = plt.cm.get_cmap('tab20')
    colors = [cmap(i) for i in range(len(labels))]
    cmap = dict(zip(labels, colors))

    # Volume of bulk myoplasm in um^3 = m^3 * 1e-18
    Vi = model.get('cell.V_i').eval()
    #Vmyo = model.get('geom.Vmyo').eval()
    #Vsl = model.get('geom.Vsl').eval()
    #Vjn = model.get('geom.Vjn').eval()
    #Vsr = model.get('geom.Vsr').eval()

    # Total internal species
    ki = (
        + d['potassium.K_i'] * (Vmyo + Vsl + Vjn)       # mM/L * L = mmol
    ) * 1e15    # mmol * 1e15 = 1e-18 mol = attomol

    nai = (
        + d['sodium.Na_i'] * Vmyo
    ) * 1e15

    cai = (
        + d['calcium.Ca_i'] * Vmyo
        + d['calcium.Ca_sl'] * Vsl
        + d['calcium.Ca_jn'] * Vjn
        + d['calcium.Ca_sr'] * Vsr
        + d['calcium.Csqn'] * Vsr
        + d['buffca.TnCL'] * Vmyo
        + d['buffca.TnCHc'] * Vmyo
        + d['buffca.TnCHm'] * Vmyo
        + d['buffca.CaM'] * Vmyo
        + d['buffca.Myoc'] * Vmyo
        + d['buffca.Myom'] * Vmyo
        + d['buffca.SRB'] * Vmyo
        + d['buffca.SLL_sl'] * Vsl
        + d['buffca.SLH_sl'] * Vsl
        + d['buffca.SLL_jn'] * Vjn
        + d['buffca.SLH_jn'] * Vjn
    ) * 1e15

    cli = (
        + d['chloride.Cl_i'] * (Vmyo + Vsl + Vjn)
    ) * 1e15


    # Compare transmembrane flux with internal amount
    time = d.time()

    # K
    flux = (dk['INaK'][-1] + dk['IK1'][-1])
    deriv = (ki[-1] - ki[-2]) / (time[-1] - time[-2])
    print(f'K flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of K at t=1000: {deriv} amol/ms')

    # Na
    flux = (dna['INaK'][-1] + dna['INaCa'][-1] + dna['INaB'][-1])
    deriv = (nai[-1] - nai[-2]) / (time[-1] - time[-2])
    print(f'Na flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of Na at t=1000: {deriv} amol/ms')

    # Ca
    flux = (dca['IpCa'][-1] + dca['INaCa'][-1] + dca['ICaB'][-1])
    deriv = (cai[-1] - cai[-2]) / (time[-1] - time[-2])
    print(f'Ca flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of Ca at t=1000: {deriv} amol/ms')

    # Cl
    flux = +(dcl['IClB'][-1])
    deriv = (cli[-1] - cli[-2]) / (time[-1] - time[-2])
    print(f'Cl flux of major currents at t=1000 : {flux} amol/ms')
    print(f'Estimated derivative of Cl at t=1000: {deriv} amol/ms')



    #
    # Create figure
    #
    print('Plotting...')
    fig = plt.figure(figsize=(18, 11))
    fig.subplots_adjust(0.048, 0.045, 0.994, 0.945)  # Left, bottom, right, top
    grid = matplotlib.gridspec.GridSpec(12, 6, wspace=0.5, hspace=0.95)

    #
    # Relative contributions per species
    #
    ax00 = fig.add_subplot(grid[0:3, 0])
    myokit.lib.plots.cumulative_current(
        dk, labels, ax00, colors=colors, integrate=False, normalise=True)
    ax00.legend(loc=(0, 1.04), ncol=10).get_frame().set_alpha(1)
    ax00.set_ylabel('Contribution to inward flux')

    ax01 = fig.add_subplot(grid[0:3, 1])
    myokit.lib.plots.cumulative_current(
        dna, labels, ax01, colors=colors, integrate=False, normalise=True)

    ax02 = fig.add_subplot(grid[0:3, 2])
    myokit.lib.plots.cumulative_current(
        dca, labels, ax02, colors=colors, integrate=False, normalise=True)

    ax03 = fig.add_subplot(grid[0:3, 3])
    myokit.lib.plots.cumulative_current(
        dcl, labels, ax03, colors=colors, integrate=False, normalise=True)

    #
    # Absolute contributions: ions / s
    #
    ax10 = fig.add_subplot(grid[3:6, 0])
    ax10.set_ylabel('Selected currents\nInward K+ flux (amol/ms)')
    sel = ['IK1', 'INaK']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dk, sel, ax10, colors=col)

    ax11 = fig.add_subplot(grid[3:6, 1])
    ax11.set_ylabel('Selected currents\nInward Na+ flux (amol/ms)')
    sel = ['INaK', 'INaB', 'INaCa']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dna, sel, ax11, colors=col)

    ax12 = fig.add_subplot(grid[3:6, 2])
    ax12.set_ylabel('Selected currents\nInward Ca2+ flux (amol/ms)')
    sel = ['IpCa', 'INaCa', 'ICaB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dca, sel, ax12, colors=col)

    ax13 = fig.add_subplot(grid[3:6, 3])
    ax13.set_ylabel('Selected currents\nInward Cl- flux (amol/ms)')
    sel = ['IClB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dcl, sel, ax13, colors=col)

    #
    # Bulk amounts
    #
    def changeplot(ax, t, y):
        ax.axhline(y[0], color='#cccccc')
        ax.plot(t, y, color='k')

    ax20 = fig.add_subplot(grid[6:8, 0])
    ax20.set_ylabel('K+ ions (amol)')
    #ax20.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
    changeplot(ax20, d.time(), ki)

    ax21 = fig.add_subplot(grid[6:8, 1])
    ax21.set_ylabel('Na+ ions (amol)')
    changeplot(ax21, d.time(), nai)

    ax22 = fig.add_subplot(grid[6:8, 2])
    ax22.set_ylabel('Ca2+ ions (amol)')
    changeplot(ax22, d.time(), cai)

    ax23 = fig.add_subplot(grid[6:8, 3])
    ax23.set_ylabel('Cl- ions (amol)')
    changeplot(ax23, d.time(), cli)

    #
    # Background current
    #
    ax31 = fig.add_subplot(grid[8:10, 1])
    ax31.set_ylabel('INaB flux (amol/s)')
    sel = ['INaB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dna, sel, ax31, colors=col)

    ax32 = fig.add_subplot(grid[8:10, 2])
    ax32.set_ylabel('ICaB flux (amol/ms)')
    sel = ['ICaB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dca, sel, ax32, colors=col)

    ax33 = fig.add_subplot(grid[8:10, 3])
    ax33.set_ylabel('[Cl]i (mM)')
    ax33.set_ylabel('IClB flux (amol/ms)')
    sel = ['IClB']
    col = [cmap[x] for x in sel]
    myokit.lib.plots.cumulative_current(dcl, sel, ax33, colors=col)

    #
    # Phase-plots for concentrations
    #
    # Set `model` to final simulation state
    model.set_state(s.state())
    model.remove_derivative_references()

    def inf(y, x):
        m = model.clone()
        x = m.get(x)
        assert(x.is_state())
        m.binding('pace').set_binding(None)
        for state in list(m.states()):
            if state != x:
                state.set_rhs(state.state_value())
                state.demote()
        return m.get(y).pyfunc()

    color1 = 'blue'
    color2 = 'brown'

    f = inf('potassium.K_i', 'potassium.K_i')
    c = model.get('potassium.K_i').state_value()
    dc = 20
    x = np.linspace(c - dc, c + dc, 100)
    ax40 = fig.add_subplot(grid[10:12, 0])
    ax40.set_xlabel('[K+]i')
    ax40.set_ylabel('dot([K+]i)')
    ax40.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
    ax40.axvline(c, color='#cccccc')
    ax40.plot(x, f(x), color=color1)

    f1 = inf('sodium.Na_jn', 'sodium.Na_jn')
    c1 = model.get('sodium.Na_jn').state_value()
    f2 = inf('sodium.Na_sl', 'sodium.Na_sl')
    c2 = model.get('sodium.Na_sl').state_value()
    c = (c1 + c2) / 2
    dc = 0.001
    x = np.linspace(c - dc, c + dc, 100)
    ax41 = fig.add_subplot(grid[10:12, 1])
    ax41.set_xlabel('[Na+]')
    ax41.set_ylabel('dot([Na+])')
    ax41.axvline(c1, color='#cccccc')
    ax41.axvline(c2, color='#cccccc', ls='--')
    ax41.axhline(0, color='#777777')
    ax41.plot(x, f1(x), color=color1, label='jn')
    ax41.plot(x, f2(x), color=color2, label='sl', ls='--')
    ax41.legend()

    f1 = inf('calcium.Ca_jn', 'calcium.Ca_jn')
    c1 = model.get('calcium.Ca_jn').state_value()
    f2 = inf('calcium.Ca_sl', 'calcium.Ca_sl')
    c2 = model.get('calcium.Ca_sl').state_value()
    c = (c1 + c2) / 2
    dc = 0.6e-4
    x = np.linspace(max(0, c - dc), c + dc, 100)
    ax42 = fig.add_subplot(grid[10:12, 2])
    ax42.set_xlabel('[Ca2+]')
    ax42.set_ylabel('dot([Ca2+])')
    ax42.axvline(c1, color='#cccccc')
    ax42.axvline(c2, color='#cccccc', ls='--')
    ax42.axhline(0, color='#777777')
    ax42.plot(x, f1(x), color=color1, label='jn')
    ax42.plot(x, f2(x), color=color2, label='sl', ls='--')

    f = inf('chloride.Cl_i', 'chloride.Cl_i')
    c = model.get('chloride.Cl_i').state_value()
    dc = 5
    x = np.linspace(c - dc, c + dc, 100)
    ax43 = fig.add_subplot(grid[10:12, 3])
    ax43.set_xlabel('[Cl-]i')
    ax43.set_ylabel('dot([Cl-]i)')
    ax43.ticklabel_format(axis='y', scilimits=(0, 0), useMathText=True)
    ax43.axvline(c, color='#cccccc')
    ax43.plot(x, f(x), color=color1)

    #
    # Influence of concentrations on currents
    #
    sub = matplotlib.gridspec.GridSpecFromSubplotSpec(
        6, 3, subplot_spec=grid[:, 4:], wspace=0.5, hspace=0.5)

    def plot_inf(ax, y, x):

        dual = False
        if x == 'membrane.V':
            ax.set_xlabel('Vr (mV)')
            dx = 4
        elif x == 'potassium.K_i':
            ax.set_xlabel('[K+]i')
            dx = 10
        elif x == 'chloride.Cl_i':
            ax.set_xlabel('[Cl-]i')
            dx = 2
        elif x == 'sodium.Na':
            dual = True
            ax.set_xlabel('[Na+] jn,sl')
            dx = 1
        elif x == 'calcium.Ca':
            dual = True
            ax.set_xlabel('[Ca2+] jn,sl')
            dx = 0.0005
        else:
            raise NotImplementedError(x)

        label = y.split('.')[1]
        ax.set_ylabel(f'{label} (pA)')
        color = cmap[y.split('.')[1]]
        #f = -1000 / 96485
        #dk[label] = d[current] * ik.get(label, 0) * f

        if dual:
            x1, y1 = x + '_jn', y + '_jn'
            x2, y2 = x + '_sl', y + '_sl'
            c1 = model.get(x1).state_value()
            c2 = model.get(x2).state_value()
            xr = max(0, min(c1 - dx, c2 - dx)), max(c1 + dx, c2 + dx)
            r = np.linspace(xr[0], xr[1], 100)
            f1 = inf(y1, x1)
            f2 = inf(y2, x2)
            ax.axvline(c1, color='#cccccc')
            ax.axvline(c2, color='#cccccc', ls='--')
            ax.axhline(f1(c1), color='#cccccc')
            ax.axhline(f2(c2), color='#cccccc', ls='--')
            ax.plot(r, f1(r), color=color, label='jn', lw=2)
            ax.plot(r, f2(r), color=color, label='sl', lw=2, ls='--')
            yc = 0.5 * (f1(c1) + f2(c2))
        else:
            c = model.get(x).state_value()
            xr = c - dx, c + dx
            r = np.linspace(xr[0], xr[1], 100)
            f = inf(y, x)
            ax.axvline(c, color='#cccccc')
            ax.axhline(f(c), color='#cccccc')
            ax.plot(r, f(r), color=color, lw=3)
            yc = f(c)

        # Set scaling
        ax.set_xlim(*xr)
        dy = 10
        ax.set_ylim(yc - dy, yc + dy)


    # IK1(K, V)
    ax = fig.add_subplot(sub[0, 0])
    plot_inf(ax, 'ik1.IK1', 'potassium.K_i')
    ax.text(1.1, 0.5, 'Positive:\n1 K+ out', transform=ax.transAxes)
    ax = fig.add_subplot(sub[0, 2])
    plot_inf(ax, 'ik1.IK1', 'membrane.V')

    # INaK(Na, V)
    ax = fig.add_subplot(sub[1, 0])
    plot_inf(ax, 'inak.INaK', 'sodium.Na')
    ax.legend(loc=(1.1, 0))
    ax.text(1.1, 0.5, 'Positive:\n3 Na+ out, 2 K+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[1, 2])
    plot_inf(ax, 'inak.INaK', 'membrane.V')

    # INaB(Na, V)
    ax = fig.add_subplot(sub[2, 0])
    plot_inf(ax, 'inab.INaB', 'sodium.Na')
    ax.legend(loc=(1.1, 0))
    ax.text(1.1, 0.5, 'Negative: 1 Na+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[2, 2])
    plot_inf(ax, 'inab.INaB', 'membrane.V')

    # INaCa(Na, Ca, V)
    ax = fig.add_subplot(sub[3, 0])
    plot_inf(ax, 'inaca.INaCa', 'sodium.Na')
    ax = fig.add_subplot(sub[3, 1])
    plot_inf(ax, 'inaca.INaCa', 'calcium.Ca')
    ax.text(-0.45, 1.07, 'Negative: 1 Ca2+ out, 3 Na+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[3, 2])
    plot_inf(ax, 'inaca.INaCa', 'membrane.V')

    # ICaB(Ca, V)
    ax = fig.add_subplot(sub[4, 0])
    plot_inf(ax, 'icab.ICaB', 'calcium.Ca')
    ax.legend(loc=(1.1, 0))
    ax.text(1.1, 0.5, 'Negative: 1 Ca2+ in', transform=ax.transAxes)
    ax = fig.add_subplot(sub[4, 2])
    plot_inf(ax, 'icab.ICaB', 'membrane.V')

    # IClB(Cl, V)
    ax = fig.add_subplot(sub[5, 0])
    plot_inf(ax, 'iclb.IClB', 'chloride.Cl_i')
    ax.text(1.1, 0.5, 'Negative: 1 Cl- out', transform=ax.transAxes)
    ax = fig.add_subplot(sub[5, 2])
    plot_inf(ax, 'iclb.IClB', 'membrane.V')

    #
    # Shading
    #
    # Add model name
    x = 0.672
    y = 0.182
    fig.patches.append(matplotlib.patches.Rectangle(
        (0, 0), 1, y, color='#eeeeee', fill=True, zorder=-1,
        transform=fig.transFigure))
    fig.patches.append(matplotlib.patches.Rectangle(
        (x, 0), 1 - x, 1, color='#eeeeee', fill=True, zorder=-1,
        transform=fig.transFigure))

    #
    # Store
    #
    path = os.path.join(fpath, fname)
    print(f'Writing to {path}')
    plt.savefig(path + '.png')
