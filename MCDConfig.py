#!/usr/bin/env python
"""Configuration file for MCDToPlanetocosmics.py, shouldn't need to be edited in normal use cases."""

#some MCD settings
dust_scenarios = {'climAVE': 1, 'climMIN': 2, 'climMAX': 3, 'dustMIN': 4, 'dustAVE': 5, 'dustMAX': 6, 'warm': 7, 'cold': 8}
MCD_DATA_DIR = '/data/etph/msl/MCD5.0/data/'
MCD_GLOBAL_CONFIG = '/data/etph/msl/MCDToPlanetocosmics/mcd_global.cfg'

#Indices of atmospheric components in extvar array
component_indices = {'CO2': 56, 'N2': 57, 'Ar': 58, 'CO': 59, 'O': 60, 'O2': 61, 'O3': 62, 'H': 63, 'H2': 64}
