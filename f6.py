#!/usr/bin/env python3
#
# Figure 6: Plots 3 INaK dependencies for the schematic
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec
import numpy as np
import os
import shared

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-6'
shared.splash(fname)

# Load model
model = shared.model('voigt')
name = model.name()
print(f'Loaded {name}')

# Prepare model
shared.prepare_model(model)

# Get some variables
inak = model.labelx('I_NaK')
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')
nai = model.labelx('Na_i')

#
# Create figure
#
fig = plt.figure(figsize=(6, 1.9))  # Two-column size, about 2/3s
fig.subplots_adjust(0.09, 0.24, 0.98, 0.98)
grid = matplotlib.gridspec.GridSpec(1, 3, wspace=0.5)


def func(model, y, x, x2=None):
    model = model.clone()
    x = model.get(x.qname())
    y = model.get(y.qname())
    for s in list(model.states()):
        if s is not x and s is not x2:
            shared.demote(s)
    if not x.is_state():
        x.promote(x.rhs().eval())
        x.set_rhs(0)
    return y.pyfunc()


# Minor ticks
matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['ytick.minor.visible'] = True

# INaK vs Ko
x = [ko]
r = (0, 15)
q = (0, 0.28)
f = func(model, inak, *x)
x = np.linspace(r[0], r[1], 100)
y = f(x)
ax1 = fig.add_subplot(grid[0, 0])
ax1.set_xlabel('[K$^+$]$_o$ (mM)')
ax1.set_ylabel('$I_{NaK}$ (A/F)')
ax1.set_xlim(*r)
ax1.set_ylim(*q)
ax1.plot(x, y)
ax1.axvline(5.4, color='#cccccc', zorder=0)

# INaK vs [Na]i
nai1 = model.get('sodium.Na_jn')
nai2 = model.get('sodium.Na_sl')
x = [nai1, nai2]
r = (0, 20)
f = func(model, inak, *x)
x = np.linspace(r[0], r[1], 100)
y = f(x)
ax2 = fig.add_subplot(grid[0, 1])
ax2.set_xlabel('[Na$^+$]$_i$ (mM)')
ax2.set_ylabel('$I_{NaK}$ (A/F)')
ax2.set_xlim(*r)
ax2.set_ylim(*q)
ax2.plot(x, y)
ax2.axvline(nai.state_value(), color='#cccccc', zorder=0)

# INaK vs Ko
x = [vm]
r = (-100, -40)
f = func(model, inak, *x)
x = np.linspace(r[0], r[1], 100)
y = f(x)
ax3 = fig.add_subplot(grid[0, 2])
ax3.set_xlabel('V (mV)')
ax3.set_ylabel('$I_{NaK}$ (A/F)')
ax3.set_xlim(*r)
ax3.set_ylim(*q)
ax3.plot(x, y)
ax3.axvline(vm.state_value(), color='#cccccc', zorder=0)


#
# Store
#
fname = f'{fname}-graphs'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.svg')
print('Done')

