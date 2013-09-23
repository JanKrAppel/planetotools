#!/usr/bin/env python

import re
from numpy import *
from matplotlib import pyplot as plt
from planetoparse import planetoparse, histdata

if not plt.isinteractive():
	plt.ion()

def __parse_title(title):
	titleparse = re.match('(.*)\s*\[(.*)\]', title)
	if titleparse is None:
		return None, title
	else:
		return titleparse.group(1), titleparse.group(2)
		
def __get_units_from_label(label):
	res = label.split(' / ')
	if len(res) == 1:
		return res[0]
	else:
		return res[1]

ENERGY_UNITS = ['MeV', 'keV', 'GeV']	
AREA_UNITS = ['cm2', 'm2', 'km2']
TIME_UNITS = ['h', 'm', 's']
COUNT_UNITS = ['nb particles']
ANGLE_UNITS = ['sr']
WEIGHT_UNITS = ['nuc']
UNIT_ORDER = ['count', 'area', 'time', 'angle', 'energy', 'weight']
UNIT_REPLACE = {'nb particles': 'particles'}

def __normalize_units(units):
	units = units.split('/')
	unitsort = {}
	unitsort['main'] = units.pop(0)
	for unit in units:
		if unit in ENERGY_UNITS:
			unitsort['energy'] = unit
		elif unit in AREA_UNITS:
			unitsort['area'] = unit
		elif unit in TIME_UNITS:
			unitsort['time'] = unit
		elif unit in COUNT_UNITS:
			unitsort['count'] = unit
		elif unit in ANGLE_UNITS:
			unitsort['angle'] = unit
		elif unit in WEIGHT_UNITS:
			unitsort['weight'] = unit
	res = unitsort['main']
	for unit in UNIT_ORDER:
		if unit in unitsort:
			res += '/' + unitsort[unit]
	for string in UNIT_REPLACE:
		res = re.sub(string, UNIT_REPLACE[string], res)
	return res
	
def __get_current_xy_labels():
	fig = plt.gcf()
	ax = fig.gca()
	return ax.get_xlabel(), ax.get_ylabel()
	
def __check_xy_units(xunits, yunits):
	match_x = __check_x_units(xunits)
	match_y = __check_y_units(yunits)
	mismatch_axes = ''
	if not match_x:
		mismatch_axes += 'X'
	if not match_y:
		if not match_x:
			mismatch_axes += ', Y'
		else:
			mismatch_axes += 'Y'
	if not (match_x or match_y):
		print 'WARNING: Units mismatch on axis', mismatch_axes
	return
	
def __check_x_units(xunits):
	cur_xunits, cur_yunits = __get_current_xy_labels()
	if not cur_xunits == '':
		cur_xunits = __get_units_from_label(cur_xunits)
		return xunits == cur_xunits
	else:
		return True

def __check_y_units(yunits):
	cur_xunits, cur_yunits = __get_current_xy_labels()
	if not cur_yunits == '':
		cur_yunits = __get_units_from_label(cur_yunits)
		return yunits == cur_yunits
	else:
		return True

class log_interpolator:
	def __init__(self, atmodata, field_x, field_y):
		self.atmodata = atmodata
		self.field_x = field_x
		self.field_y = field_y
		return
		
	def __call__(self, x):
		return 10**interp(x, self.atmodata.data[self.field_x], log10(self.atmodata.data[self.field_y]))

def convert_edep_to_LET(profile, atmodata):
	"""Convert a energy deposition profile in rad/s vs km to keV/um vs km. Options needed are the deposition profile to convert and the atmodata/mcddata atmosphere data used in the simulation run."""
	res = histdata(copyhist = profile)
	tmp, yunits = __parse_title(profile.params['Xaxis'])
	if yunits == 'km':
		x_field = 'xz'
	elif yunits == 'g/cm2':
		x_field = 'shield_depth'
	else:
		print 'ERROR: Unknown altitude scaling, aborting'
		return None
	interpolate = log_interpolator(atmodata, x_field, 'dens')
	for line in res.data:
		dens = interpolate(line[2])
		line[3:] /= 1e3 #in J*cm**2/g
		#convert to J/cm:
		line[3:] *= dens
		#convert to keV/cm:
		line[3:] *= 6.241509e15
		#convert to keV/um:
		line[3:] /= 1e4
	title, xunits = __parse_title(res.params['Title'])
	title = re.sub('Deposited energy', 'LET', title)
	res.params['Title'] = title + ' [keV/um]'
	return res
		
def plot_edep_profile(hist, *args, **kwargs):
	"""Plots energy deposition profiles. Pass the profile as available through planetoparse to plot, additional arguments are passed to the Matplotlib plotting function (errorbar)."""
	if hist.isempty():
		print 'WARNING: Unable to plot, histogram is all-zero.'
		return
	if not 'capsize' in kwargs:
		capsize = 0
	else:
		capsize = kwargs['capsize']
		kwargs.pop('capsize')
	bin_width = hist.data[:,1] - hist.data[:,0]
	plt.errorbar(hist.data[:,3] / bin_width, hist.data[:,2], xerr = hist.data[:,4] / bin_width, marker='.', capsize = capsize, *args, **kwargs)
	title, xunits = __parse_title(hist.params['Title'])
	ylabel, yunits = __parse_title(hist.params['Xaxis'])
	__check_xy_units(xunits, yunits)
	plt.title(title)
	plt.xlabel('Deposited energy / ' + xunits)
	plt.ylabel(ylabel + ' / ' + yunits)
	plt.xscale('log')
	plt.ylim(amin(hist.data[:,2]), amax(hist.data[:,2]))
	plt.show(block = False)
	return
	
def plot_1d_hist(hist, scale_by = 1., label_detector = False, scale_by_width = True, *args, **kwargs):
	"""Plots 1D histograms. Pass the histogram as available through planetoparse to plot, additional arguments are passed to the Matplotlib plotting function (errorbar)."""
	if hist.isempty():
		print 'WARNING: Unable to plot, histogram is all-zero.'
		return
	if not 'label' in kwargs:
		label = hist.particle
		if label_detector:
			label += ', detector ' + str(hist.detector)
	else:
		label = kwargs['label']
		kwargs.pop('label')
	if not 'capsize' in kwargs:
		capsize = 0
	else:
		capsize = kwargs['capsize']
		kwargs.pop('capsize')
	if scale_by_width:
		bin_width = hist.data[:,1] - hist.data[:,0]
		if (bin_width == 0.).all():
			print 'WARNING: Unable to scale by bin width'
			scale_by_width = False
			bin_width = ones(len(bin_width))
	else:
		bin_width = ones(len(hist.data))
	plt.errorbar(hist.data[:,2], hist.data[:,3] * scale_by / bin_width, xerr = bin_width / 2, yerr = hist.data[:,4] * scale_by / bin_width, marker='.', label = label, capsize = capsize, *args, **kwargs)
	title, units = __parse_title(hist.params['Title'])
	plt.title(title)
	xlabel, xunits = __parse_title(hist.params['Xaxis'])
	xunits = __normalize_units(xunits)
	if scale_by_width:
		yunits = __normalize_units(units + '/' + xunits)
	else:
		yunits = __normalize_units(units)
	__check_xy_units(xunits, yunits)
	plt.xlabel(xlabel + ' / ' + xunits)
	plt.ylabel(yunits)
	plt.xscale('log')
	plt.yscale('log')
	plt.legend(loc = 'best')
	plt.show(block = False)
	return
	
def plot_array_hist(array, scale_by = 1., *args, **kwargs):
	"""Plots 1D histograms from numpy arrays. Additional arguments are passed to the Matplotlib plotting function (errorbar)."""
	if not 'capsize' in kwargs:
		capsize = 0
	else:
		capsize = kwargs['capsize']
		kwargs.pop('capsize')
	#errors are included:
	if array.shape[1] == 5:
		scale_by_width = True
		bin_width = array[:,1] - array[:,0]
		if (bin_width == 0.).all():
			print 'WARNING: Unable to scale by bin width'
			scale_by_width = False
			bin_width = ones(len(bin_width))
		plt.errorbar(array[:,2], array[:,3] * scale_by / bin_width, xerr = bin_width / 2, yerr = array[:,4] * scale_by / bin_width, marker='.', capsize = capsize, *args, **kwargs)
	#errors are not included, neither are binwidths:
	elif array.shape[1] == 2:
		plt.plot(array[:,0], array[:,1], *args, **kwargs)
	plt.xscale('log')
	plt.yscale('log')
	plt.show(block = False)
	return
	
def scale_array_per_nuc(array, weight):
	"""Scales an array that contains Planetocosmics result data per nucleus. Requires the following columns:
	1. lower bin edge 2. upper bin edge 3. bin middle 4. bin height 5. bin error"""
	array[:,:3] /= weight
	array[:,3:] *= weight
	return

def plot_2d_hist(hist, logscale = True, *args, **kwargs):
	"""Plots 2D histogram data. Pass the histogram as available through planetoparse to plot."""
	if hist.isempty():
		print 'WARNING: Unable to plot, histogram is all-zero.'
		return
	histdat = []
	xedges = []
	yedges = []
	for line in hist.data:
		if not line[0] in xedges:
			xedges.append(line[0])
		if not line[1] in xedges:
			xedges.append(line[1])
		if not line[2] in yedges:
			yedges.append(line[2])
		if not line[3] in yedges:
			yedges.append(line[3])
		if logscale:
			histdat.append(log10(line[4]))
		else:
			histdat.append(line[4])
	histdat = array(histdat).reshape((len(xedges) - 1, len(yedges) - 1))
	histdat[histdat == -inf] = 0.
	xedges = array(xedges)
	yedges = array(yedges)
	xedges.sort()
	yedges.sort()
	plt.pcolormesh(xedges, yedges, histdat.transpose(), *args, **kwargs)
	plt.xlabel(hist.params['Xaxis'])
	plt.ylabel(hist.params['Yaxis'])
	fig = plt.gcf()
	ax = fig.get_axes()
	title, units = __parse_title(hist.params['Title'])
	#remove colorbar, if already present
	if len(ax) == 2:
		if ax[0].get_title() == title and ax[1].get_title() == '':
			fig.delaxes(ax[1])
			fig.subplots_adjust()
	cbar=plt.colorbar()
	plt.title(title)
	if logscale:
		units = 'log ' + units
	cbar.set_label(units)
	plt.show(block = False)
	return
	
def plot_cosmonuc(results, logscale = True, *args, **kwargs):
	"""Plots 2D histogram of cosmogenic nuclides. Pass the planetoparse instance to plot."""
	if results.cosmonuc.isempty():
		print 'WARNING: Unable to plot, histogram is all-zero.'
		return
	if not 's' in kwargs:
		s = 150
	else:
		s = kwargs['s']
		kwargs.pop('s')
	if not 'marker' in kwargs:
		marker = 's'
	else:
		marker = kwargs['marker']
		kwargs.pop('marker')
	if not 'lw' in kwargs:
		lw = 0
	else:
		lw = kwargs['lw']
		kwargs.pop('lw')
	if logscale:
		c = log10(results.cosmonuc.data[:,4])
	else:
		c = results.cosmonuc.data[:,4]
	plt.scatter(results.cosmonuc.data[:,0] + .5, results.cosmonuc.data[:,2] + .5, c = c, s = s, marker = marker, lw = lw, *args, **kwargs)
	plt.xlabel(results.cosmonuc.params['Xaxis'])
	plt.ylabel(results.cosmonuc.params['Yaxis'])
	plt.xlim([0,25])
	plt.xticks(arange(0, 26, 2))
	plt.ylim([0,25])
	plt.yticks(arange(0, 26, 2))
	fig = plt.gcf()
	ax = fig.get_axes()
	title, units = __parse_title(results.cosmonuc.params['Title'])
	#remove colorbar, if already present
	if len(ax) == 2:
		if ax[0].get_title() == title and ax[1].get_title() == '':
			fig.delaxes(ax[1])
			fig.subplots_adjust()
	cbar=plt.colorbar()
	plt.title(title)
	if logscale:
		units = 'log ' + units
	cbar.set_label(units)
	plt.show(block = False)
	return
	
def plot_detector_levels(fluxhists, plot_only = [], dont_plot = []):
	"""Plots the detector levels of the given flux histograms into the current plot."""
	left, right = plt.xlim()
	if not plot_only == []:
		detectors = plot_only
	else:
		detectors = fluxhists.keys()
	for detector in dont_plot:
		if detector in detectors:
			detectors.remove(detector)
	altitudes = []
	for detector in detectors:
		altitude, alt_unit = fluxhists[detector].params['Altitude'].split(' ')
		altitude = float64(altitude)
		altitudes.append(altitude)
		plt.plot((left, right), (altitude, altitude), 'k')
		plt.text(left, altitude + 1., 'Detector ' + str(detector))
	if not __check_y_units(alt_unit):
		print 'WARNING: Units mismatch on axis Y'
	plt.ylim(amin(altitudes), amax(altitudes) + 10)
	plt.ylabel('Altitude / ' + alt_unit)
	plt.show(block = False)
	return
	
def combine_histograms(*args):
	"""Combines the different histograms for multiple planetoparse instances.
	WARNING: Only the primaries count, cosmonuc 2D histograms, primary particle fluxes and up/down flux histograms are combined!"""
	if len(args) > 1:
		res = planetoparse()
		for addthis in args:
			#check and set normalisation
			if not res.normalisation == '':
				if not res.normalisation == addthis.normalisation:
					print 'WARNING: Skipping result, different normalisation methods'
					print '\t' + addthis.normalisation + ' vs. ' + res.normalisation + ' as expected'
					continue
			else:
				res.normalisation += addthis.normalisation
			#combine primary counts
			res.primaries += addthis.primaries
			#combine number of events counts
			if 'nb_of_events' in res.params:
				res.params['nb_of_events'] += addthis.params['nb_of_events']
			else:
				res.params['nb_of_events'] = float64(0.) + addthis.params['nb_of_events']
			#combine cosmogenic nuclides
			if not res.cosmonuc is None:
				res.cosmonuc.data[:,4] += addthis.cosmonuc.data[:,4]
			else:
				res.cosmonuc = histdata(copyhist = addthis.cosmonuc)
			#combine upward fluxes
			for particle in addthis.flux_up:
				for detector in addthis.flux_up[particle]:
					if not particle in res.flux_up:
						res.flux_up[particle] = {}
					if detector in res.flux_up[particle]:
						res.flux_up[particle][detector] = __combine_single_hists(res.flux_up[particle][detector], addthis.flux_up[particle][detector])
					else:
						res.flux_up[particle][detector] = histdata(copyhist = addthis.flux_up[particle][detector])
			#combine downward fluxes
			for particle in addthis.flux_down:
				for detector in addthis.flux_down[particle]:
					if not particle in res.flux_down:
						res.flux_down[particle] = {}
					if detector in res.flux_down[particle]:
						res.flux_down[particle][detector] = __combine_single_hists(res.flux_down[particle][detector], addthis.flux_down[particle][detector])
					else:
						res.flux_down[particle][detector] = histdata(copyhist = addthis.flux_down[particle][detector])
			#combine primary fluxes
			for particle in addthis.primhists:
				if not particle in res.primhists:
					res.primhists[particle] = histdata(copyhist = addthis.primhists[particle])
				else:
					res.primhists[particle] = __combine_single_hists(res.primhists[particle], addthis.primhists[particle])
		return res
	else:
		return args[0]
		
def __combine_single_hists(hist1, hist2):
	from planetoparse import histdata
	res = histdata(copyhist = hist1)
	if not (res.data[:,0] == hist2.data[:,0]).all() or not (res.data[:,1] == hist2.data[:,1]).all() or not (res.data[:,2] == hist2.data[:,2]).all():
		print 'ERROR: Unable to combine histograms, binning is different.'
		return hist1
	xunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
	yunits1 = __normalize_units(__parse_title(res.params['Title'])[1])
	xunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
	yunits2 = __normalize_units(__parse_title(hist2.params['Title'])[1])
	if not (xunits1 == xunits2 and yunits1 == yunits2):
		print 'ERROR: Unable to combine histograms, units mismatch.'
	res.data[:,3] += hist2.data[:,3]
	res.data[:,4] = sqrt(res.data [:,4]**2 + hist2.data[:,4]**2)
	return res
		
def plot_primaries(results, *args, **kwargs):
	"""Plot all primary particle fluxes in a result."""
	for particle in results.primhists:
		if particle == 'alpha':
			results.primhists[particle].scale_per_nuc(4)
		plot_1d_hist(results.primhists[particle], *args, **kwargs)
	return
	
def plot_proton_alpha(results, detector, scale_per_nuc = True, *args, **kwargs):
	"""Plot downward proton and alpha fluxes at the given detector."""
	plot_1d_hist(results.flux_down['proton'][detector], *args, **kwargs)
	if scale_per_nuc:
		results.flux_down['alpha'][detector].scale_per_nuc(4)
	plot_1d_hist(results.flux_down['alpha'][detector], *args, **kwargs)
	return
	
def plot_neutrals(results, detector, *args, **kwargs):
	"""Plot downward neutron and gamma fluxes at the given detector."""
	plot_1d_hist(results.flux_down['neutron'][detector], *args, **kwargs)
	plot_1d_hist(results.flux_down['gamma'][detector], *args, **kwargs)
	return
	
def plot_proton_alpha_comparison(results, detector, scale_per_nuc = True, *args, **kwargs):
	"""Plot primary, proton/alpha and neutral downward fluxes at detector."""
	plot_primaries(results, *args, **kwargs)
	plot_proton_alpha(results, detector, scale_per_nuc = scale_per_nuc, *args, **kwargs)
	return
	
def plot_cosmonuc_comparison(results1, results2, label1 = None, label2 = None, legend = True, *args, **kwargs):
	"""Plot a comparison of two cosmogenic nuclide histograms."""
	for entry in ['label', 'lw', 's', 'marker', 'edgecolor', 'legend']:
		if entry in kwargs:
			kwargs.pop(entry)
	if label1 is None:
		plot_cosmonuc(results1, *args, **kwargs)
	else:
		plot_cosmonuc(results1, label = label1, *args, **kwargs)
	if label2 is None:
		plot_cosmonuc(results2, marker = 'o', lw = 1, edgecolor = 'w', s = 100, *args, **kwargs)
	else:
		plot_cosmonuc(results2, label = label2, marker = 'o', lw = 1, edgecolor = 'w', s = 100, *args, **kwargs)
	if legend:
		plt.legend(loc = 'upper left')
	plt.title('Comparison of cosmogenic nuclide production')
	return
	
def get_normalization_factors(results):
	"""Get a dictionary with the normalization factors for all particles in the flux_up and flux_down sections of the given results."""
	res = {}
	for flux_list in [results.flux_down, results.flux_up]:
		for particle in flux_list:
			detector = flux_list[particle].keys()[-1] #this is a completely arbitrary choice
			if 'normalisation_factor' in flux_list[particle][detector].params:
				factor = flux_list[particle][detector].params['normalisation_factor']
			else:
				factor = nan
			if not particle in res:
				res[particle] = factor
	return res
	
def print_list_info(histlist, indices = None):
	"""Prints information on histograms in the given histogram list. Optional second argument selects the indices to print."""
	if indices is None or type(indices) == dict:
		indices = arange(0, len(histlist))
	if not hasattr(histlist, '__iter__'):
		try:
			indices = int(indices)
		except ValueError:
			print 'WARNING: Unable to convert given index, defaulting to all' 
			indices = arange(0, len(histlist))
	message = ''
	for index in indices:
		message += '{index: >3d}'.format(index = index) + ': {particle: <10s}'.format(particle = histlist[index].particle) + ', detector {detector: >3d}'.format(detector = histlist[index].detector) + ': '
		message += histlist[index].params['Title'] + '\n'
	print message[:-1]
	return

def integrate_fluxhist(histogram, limits = None):
	"""Integrates the flux stored in the given 1D histogram. Pass a tuple (Emin, Emax) as second argument to set the integration range."""
	if isinstance(histogram, histdata):
		if not limits is None:
			mask = (histogram.data[:,0] >= limits[0]) * (histogram.data[:,1] <= limits[1])
		else:
			mask = ones(len(histogram.data), dtype = bool)
		return sum(histogram.data[:,3][mask])
	elif type(histogram) == dict:
		if not limits is None:
			mask = (histogram['x'] >= limits[0]) * (histogram['x'] <= limits[1])
		else:
			mask = ones(len(histogram['x']), dtype = bool)
		return sum(histogram['y'][mask] * histogram['bin_widths'][mask])
		
def integrate_2d_fluxhist(histogram, xlimits = None, ylimits = None):
	"""Integrates the flux stored in the given 2D histogram. Pass a tuple (Emin, Emax) as xlimits or ylimits keyword arguments to set the integration range."""
	if not xlimits is None:
		xmask = (histogram.data[:,0] >= xlimits[0]) * (histogram.data[:,1] <= xlimits[1])
	else:
		xmask = ones(len(histogram.data), dtype = bool)
	if not ylimits is None:
		ymask = (histogram.data[:,2] >= ylimits[0]) * (histogram.data[:,3] <= ylimits[1])
	else:
		ymask = ones(len(histogram.data), dtype = bool)
	mask = xmask * ymask
	return sum(histogram.data[:,4][mask])

