#!/usr/bin/env python

import re
from numpy import *
from matplotlib import pyplot as plt
from planetoparse import planetoparse

def __parse_title(title):
	titleparse = re.match('(.*)\s*\[(.*)\]', title)
	return titleparse.group(1), titleparse.group(2)	

def plot_1d_hist(hist, label_detector = False, *args, **kwargs):
	"""Plots 1D histograms. Pass the histogram as available through planetoparse to plot, additional arguments are passed to the Matplotlib plotting function (errorbar)."""
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
	bin_width = hist.data[:,1] - hist.data[:,0]
	plt.errorbar(hist.data[:,2], hist.data[:,3] / bin_width, xerr = bin_width / 2, yerr = hist.data[:,4] / bin_width, marker='.', label = label, capsize = capsize, *args, **kwargs)
	title, units = __parse_title(hist.params['Title'])
	plt.title(title)
	plt.xlabel(hist.params['Xaxis'])
	tmp, yunit = __parse_title(hist.params['Xaxis'])
	plt.ylabel(units + '/' + yunit)
	plt.xscale('log')
	plt.yscale('log')
	plt.legend()
	plt.show()
	return
	
def plot_2d_hist(hist, *args, **kwargs):
	"""Plots 2D histogram of cosmogenic nuclides. Pass the histogram as available through planetoparse to plot."""
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
	plt.scatter(hist.data[:,0] + .5, hist.data[:,2] + .5, c = log10(hist.data[:,4]), s = 150, marker = 's', lw = 0)
	plt.xlabel(hist.params['Xaxis'])
	plt.ylabel(hist.params['Yaxis'])
	plt.xlim([0,25])
	plt.xticks(arange(0, 26, 2))
	plt.ylim([0,25])
	plt.yticks(arange(0, 26, 2))
	cbar=plt.colorbar()
	title, units = __parse_title(hist.params['Title'])
	plt.title(title)
	cbar.set_label('log ' + units)
	plt.show()
	return
	
def combine_histograms(*args):
	"""Combines the different histograms for multiple planetoparse instances.
	WARNING: Only the primaries count, cosmonuc 2D histograms, and up/down flux histograms are combined!"""
	if len(args) > 1:
		res = planetoparse()
		for addthis in args:
			res.primaries += addthis.primaries
			#combine cosmogenic nuclides
			if not res.cosmonuc is None:
				res.cosmonuc.data[:,4] += addthis.cosmonuc.data[:,4]
			else:
				res.cosmonuc = addthis.cosmonuc
			#combine upward fluxes
			for particle in addthis.flux_up:
				for detector in addthis.flux_up[particle]:
					if not particle in res.flux_up:
						res.flux_up[particle] = {}
					if detector in res.flux_up[particle]:
						res.flux_up[particle][detector] = __combine_single_hists(res.flux_up[particle][detector], addthis.flux_up[particle][detector])
					else:
						res.flux_up[particle][detector] = addthis.flux_up[particle][detector]
			#combine downward fluxes
			for particle in addthis.flux_down:
				for detector in addthis.flux_down[particle]:
					if not particle in res.flux_down:
						res.flux_down[particle] = {}
					if detector in res.flux_down[particle]:
						res.flux_down[particle][detector] = __combine_single_hists(res.flux_down[particle][detector], addthis.flux_down[particle][detector])
					else:
						res.flux_down[particle][detector] = addthis.flux_down[particle][detector]
			#combine primary fluxes
			for particle in addthis.primhists:
				if not particle in res.primhists:
					res.primhists[particle] = addthis.primhists[particle]
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
	else:
		res.data[:,3] += hist2.data[:,3]
		res.data[:,4] += hist2.data[:,4]
		return res
		
def plot_primaries(results, *args, **kwargs):
	"""Plot all primary particle fluxes in a result."""
	for particle in results.primhists:
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
			
if __name__ == '__main__':
	pass
