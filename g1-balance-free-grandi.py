#!/usr/bin/env python3
#
# Ionic balance in the freshly unclamped Grandi model
#
import shared
import shared_balance

fpath = shared.figure_dir()
fname = 'figure-g1-balance-free-grandi'
model = 'free-grandi-2011.mmt'

shared_balance.grandi(fpath, fname, model)
