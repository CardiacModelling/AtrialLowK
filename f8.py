#!/usr/bin/env python3
#
# Plot 3 INaK dependencies for the schematic in figure 8
#
#import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.gridspec
import myokit
#import numpy as np
#import os
import shared

# Get path for figures
fpath = shared.figure_dir()
fname = 'figure-8'

# Load model
model = shared.model('voigt')
name = model.name()
print(f'Loaded {name}')

# Prepare model
shared.prepare_model(model, fix_cleft_ko=True, pre_pace=False)

# Get some variables
inak = model.labelx('I_NaK')
vm = model.labelx('membrane_potential')
ko = model.labelx('K_o')
nai = model.labelx('Na_i')

#
# Create figure
#
fig = plt.figure(figsize=(9, 4))  # Two-column size
fig.subplots_adjust(0.10, 0.060, 0.98, 0.92)
grid = matplotlib.gridspec.GridSpec(3, 1, wspace=0.35)

# Add model name
'''
name_font = {
    'fontsize': 12,
    'verticalalignment': 'center',
    'horizontalalignment': 'center',
}
fig.text(0.5, 0.98, shared.fancy_name(model), name_font)

# Add panel letters
letter_font = {'weight': 'bold', 'fontsize': 16}
fig.text(0.003, 0.940, 'A', letter_font)
fig.text(0.003, 0.620, 'B', letter_font)
fig.text(0.003, 0.500, 'C', letter_font)
fig.text(0.003, 0.380, 'D', letter_font)
fig.text(0.003, 0.260, 'E', letter_font)
fig.text(0.003, 0.140, 'F', letter_font)
'''


#
# Store
#
fname = f'{fname}-graphs'
path = os.path.join(fpath, fname)
print(f'Writing figure to {path}')
plt.savefig(f'{path}.pdf')

