#!/usr/bin/env python

from numpy import *
from matplotlib import pyplot as plt
from mcdparse import *
from planetotools import __check_xy_units

if not plt.isinteractive():
	plt.ion()

def plot_temp_profile(mcddata, *args, **kwargs):
	"""Plot a temperature profile for the given MCD data object."""
	plt.plot(mcddata.data['temp'], mcddata.data['xz'], *args, **kwargs)
	__check_xy_units('K', 'km')
	plt.xlabel('Temperature / K')
	plt.ylabel('Height / km')
	plt.ylim(mcddata.params['surface_height'], mcddata.params['max_height'])
	plt.title('Temperature profile, ' + mcddata.get_coords_scenario_str())
	plt.show(block = False)
	return
	
def plot_pres_profile(mcddata, *args, **kwargs):
	"""Plot a pressure profile for the given MCD data object."""
	plt.plot(mcddata.data['pres'], mcddata.data['xz'], *args, **kwargs)
	__check_xy_units('Pa', 'km')
	plt.xlabel('Pressure / Pa')
	plt.ylabel('Height / km')
	plt.xscale('log')
	plt.ylim(mcddata.params['surface_height'], mcddata.params['max_height'])
	plt.title('Pressure profile, ' + mcddata.get_coords_scenario_str())
	plt.show(block = False)
	return
		
def plot_dens_profile(mcddata, *args, **kwargs):
	"""Plot a density profile for the given MCD data object."""
	plt.plot(mcddata.data['dens'], mcddata.data['xz'], *args, **kwargs)
	__check_xy_units('g/cm3', 'km')
	plt.xlabel('Density / g/cm3')
	plt.ylabel('Height / km')
	plt.xscale('log')
	plt.ylim(mcddata.params['surface_height'], mcddata.params['max_height'])
	plt.title('Density profile, ' + mcddata.get_coords_scenario_str())
	plt.show(block = False)
	return

def plot_comp_profile(mcddata, plot_only = None, dont_plot = None, *args, **kwargs):
	"""Plot a composition profile for the given MCD data object."""
	if 'label' in kwargs:
		kwargs.pop('label')
	if not plot_only is None:
		components = plot_only
	else:
		components = mcddata.data['comp'].keys()
	if not dont_plot is None:
		for component in dont_plot:
			if component in components:
				components.remove(component)
	for component in components:
		plt.plot(mcddata.data['comp'][component], mcddata.data['xz'], label = component, *args, **kwargs)
	plt.legend(loc = 'best')
	__check_xy_units('Volume mixing ratio', 'km')
	plt.xlabel('Volume mixing ratio')
	plt.ylabel('Height / km')
	plt.xscale('log')
	plt.ylim(mcddata.params['surface_height'], mcddata.params['max_height'])
	plt.title('Composition profile, ' + mcddata.get_coords_scenario_str())
	plt.show(block = False)
	return
	
def plot_shield_depth_profile(mcddata, *args, **kwargs):
	"""Plot a shielding depth profile for the given MCD data object."""
	plt.plot(mcddata.data['shield_depth'], mcddata.data['xz'], *args, **kwargs)
	plt.xlabel('Shielding depth / g/cm2')
	plt.ylabel('Height / km')
	plt.xscale('log')
	plt.ylim(mcddata.params['surface_height'], mcddata.params['max_height'])
	plt.title('Shielding depth profile, ' + mcddata.get_coords_scenario_str())
	plt.show(block = False)
	return
	
def plot_data_overview(mcddata, title = False, *args, **kwargs):
	"""Plot a data overview for the given MCD data object.
	This creates a subplot with temperature, pressure, density and composition profiles."""
	plt.subplot(221)
	plot_temp_profile(mcddata, *args, **kwargs)
	plt.subplot(222)
	plot_pres_profile(mcddata, *args, **kwargs)
	plt.subplot(223)
	plot_dens_profile(mcddata, *args, **kwargs)
	plt.subplot(224)
	plot_shield_depth_profile(mcddata, *args, **kwargs)
	if title:
		plt.suptitle('Data overview for ' + mcddata.params['filename'])
	pass

