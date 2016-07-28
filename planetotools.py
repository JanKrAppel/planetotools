#!/usr/bin/env python
"""Contains plotting and data analysis tools for working with planetoparse
instances."""

import re
from numpy import *
plotting_available = True
try:
    from matplotlib import pyplot as plt
except ImportError:
    print 'WARNING: matplotlib not available, plotting functions not available'
    plotting_available = False
from planetoparse import *

if plotting_available:
    if not plt.isinteractive():
        plt.ion()

def __parse_title(title):
    """Parses plot title and unit information from planetoparse plot titles.
    """
    titleparse = re.match('(.*)\s*\[(.*)\]', title)
    if titleparse is None:
        return None, title
    else:
        return titleparse.group(1), titleparse.group(2)
        
def __get_units_from_label(label):
    """Returns unit informations from plot labels."""
    res = label.split(' / ')
    if len(res) == 1:
        return res[0]
    else:
        return res[1]

ENERGY_UNITS = ['MeV', 'keV', 'GeV', 'eV']    
AREA_UNITS = ['cm2', 'm2', 'km2']
TIME_UNITS = ['h', 'm', 's']
COUNT_UNITS = ['nb particles', '#']
ANGLE_UNITS = ['sr']
WEIGHT_UNITS = ['nuc']
UNIT_ORDER = ['count', 'area', 'time', 'angle', 'energy', 'weight']
UNIT_REPLACE = {'nb particles': 'particles'}
ENERGY_SCALERS = {None: 1, 'GeV': 0.001, 'MeV': 1, 'eV': 1000000.0,
                  'keV': 1000.0}

def __normalize_units(units):
    """Returns unit strings in a normalized order."""
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
    """Returns X and Y axis labels of the current plot figure."""
    fig = plt.gcf()
    ax = fig.gca()
    return ax.get_xlabel(), ax.get_ylabel()
    
def __check_xy_units(xunits, yunits):
    """Checks if X and Y units match the current plot."""
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
    """Checks if X units match the current plot."""
    cur_xunits, cur_yunits = __get_current_xy_labels()
    if not cur_xunits == '':
        cur_xunits = __get_units_from_label(cur_xunits)
        return xunits == cur_xunits
    else:
        return True

def __check_y_units(yunits):
    """Checks if Y units match the current plot."""
    cur_xunits, cur_yunits = __get_current_xy_labels()
    if not cur_yunits == '':
        cur_yunits = __get_units_from_label(cur_yunits)
        return yunits == cur_yunits
    else:
        return True

class log_interpolator:
    """Interpolates data on a log scale."""
    def __init__(self, atmodata, field_x, field_y):
        self.atmodata = atmodata
        self.field_x = field_x
        self.field_y = field_y
        return
        
    def __call__(self, x):
        return 10**interp(x, self.atmodata.data[self.field_x], 
                          log10(self.atmodata.data[self.field_y]))

if plotting_available:
    def plot_edep_profile(hist, errorbars=True, *args, **kwargs):
        """Plots energy deposition profiles. Pass the profile as available 
        through planetoparse to plot, additional arguments are passed to the 
        Matplotlib plotting function (errorbar)."""
        if hist.isempty():
            print 'WARNING: Unable to plot, histogram is all-zero.'
            return
        if not 'capsize' in kwargs:
            capsize = 0
        else:
            capsize = kwargs['capsize']
            kwargs.pop('capsize')
        bin_width = hist.data[:, 1] - hist.data[:, 0]
        if errorbars:
            plt.errorbar(hist.data[:, 3] / bin_width, hist.data[:, 2], 
                         xerr = hist.data[:,4] / bin_width, marker='.', 
                         capsize = capsize, *args, **kwargs)
        else:
            plt.plot(hist.data[:, 3] / bin_width, hist.data[:, 2], marker='.', 
                     *args, **kwargs)
        title, xunits = __parse_title(hist.params['Title'])
        ylabel, yunits = __parse_title(hist.params['Xaxis'])
        __check_xy_units(xunits, yunits)
        plt.title(title)
        plt.xlabel('Deposited energy / ' + xunits)
        plt.ylabel(ylabel + ' / ' + yunits)
        plt.xscale('log')
        plt.ylim(amin(hist.data[:, 2]), amax(hist.data[:, 2]))
        plt.show(block = False)
        return
    
if plotting_available:
    def plot_1d_hist(hist, scale_by = 1., label_detector = False, 
                     scale_by_width = True, xlims = (-inf, inf), 
                     energy_scale = None, errorbars = True, *args, **kwargs):
        """Plots 1D histograms. Pass the histogram as available through 
        planetoparse to plot, additional arguments are passed to the Matplotlib 
        plotting function (errorbar)."""
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
        mask = (hist.data[:, 2] > xlims[0]) * (hist.data[:, 2] < xlims[1])
        data = hist.data[mask]
        if not energy_scale in ENERGY_SCALERS:
            print 'WARNING: Invalid energy scale specified, defaulting to MeV'
            energy_scale = None
        data[:, :3] *= ENERGY_SCALERS[energy_scale]
        params = {}
        for entry in hist.params:
            params[entry] = hist.params[entry]
        if not energy_scale is None:
            tmp = __parse_title(params['Xaxis'])
            if not tmp[0] is None:
                params['Xaxis'] = tmp[0] + '[' + energy_scale + ']'
        if scale_by_width:
            bin_width = data[:, 1] - data[:, 0]
            if (bin_width == 0.).all():
                print 'WARNING: Unable to scale by bin width'
                scale_by_width = False
                bin_width = ones(len(bin_width))
        else:
            bin_width = ones(len(data))
        if errorbars:
            plt.errorbar(data[:, 2], data[:, 3] * scale_by / bin_width,
                         xerr = bin_width / 2, 
                         yerr = data[:, 4] * scale_by / bin_width, marker='.',
                         label = label, capsize = capsize, *args, **kwargs)
        else:
            plt.plot(data[:, 2], data[:, 3] * scale_by / bin_width,
                     label = label, *args, **kwargs)            
        title, units = __parse_title(params['Title'])
        plt.title(title)
        xlabel, xunits = __parse_title(params['Xaxis'])
        if xlabel is None:
            xlabel = params['Xaxis']
        else:
            xunits = __normalize_units(xunits)
            xlabel += ' / ' + xunits
        if scale_by_width and not xunits is None:
            yunits = __normalize_units(units + '/' + xunits)
        elif scale_by_width and xunits is None:
            yunits = __normalize_units(units)
        else:
            yunits = __normalize_units(units)
        __check_xy_units(xunits, yunits)
        plt.xlabel(xlabel)
        plt.ylabel(yunits)
        if __is_logscale(hist.data[:, 2]):
            plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc = 'best')
        plt.show(block = False)
        return
    
if plotting_available:
    def plot_array_hist(array, scale_by = 1., *args, **kwargs):
        """Plots 1D histograms from numpy arrays. Additional arguments are passed
        to the Matplotlib plotting function (errorbar)."""
        if not 'capsize' in kwargs:
            capsize = 0
        else:
            capsize = kwargs['capsize']
            kwargs.pop('capsize')
        #errors are included:
        if array.shape[1] == 5:
            bin_width = array[:, 1] - array[:, 0]
            if (bin_width == 0.).all():
                print 'WARNING: Unable to scale by bin width'
                bin_width = ones(len(bin_width))
            plt.errorbar(array[:, 2], array[:, 3] * scale_by / bin_width, 
                         xerr = bin_width / 2, 
                         yerr = array[:, 4] * scale_by / bin_width, 
                         marker='.', capsize = capsize, *args, **kwargs)
        #errors are not included, neither are binwidths:
        elif array.shape[1] == 2:
            plt.plot(array[:, 0], array[:, 1], *args, **kwargs)
        if __is_logscale(array[:, 2]):
            plt.xscale('log')
        plt.yscale('log')
        plt.show(block = False)
        return
    
def scale_array_per_nuc(array, weight):
    """Scales an array that contains Planetocosmics result data per nucleus.
    Requires the following columns:
    1. lower bin edge
    2. upper bin edge
    3. bin middle
    4. bin height
    5. bin error"""
    array[:, :3] /= weight
    array[:, 3:] *= weight
    return

if plotting_available:
    def plot_2d_hist(hist, logscale = True, scale_by_widths = False,
                     *args, **kwargs):
        """Plots 2D histogram data. Pass the histogram as available through 
        planetoparse to plot."""
        if hist.isempty():
            print 'WARNING: Unable to plot, histogram is all-zero.'
            return
        if scale_by_widths:
            #scale flux by bin width /MeV
            data = empty_like(hist.data)
            data[:] = hist.data
            data[:, 4] /= data[:, 1] - data[:, 0]
        else:
            data = hist.data
        histdat = []
        xedges = []
        yedges = []
        for line in data:
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
        if __is_logscale(xedges):
            plt.xscale('log')
        if __is_logscale(yedges):
            plt.yscale('log')
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
        cbar = plt.colorbar()
        plt.title(title)
        if scale_by_widths:
            units += '/' + __parse_title(hist.params['Xaxis'])[1]
        if logscale:
            units = 'log ' + units
        cbar.set_label(__normalize_units(units))
        plt.show(block = False)
        return
    
def __is_logscale(axis):
    """Checks whether or not a given axis is logscale."""
    if (diff(diff(axis)) <= 1e-15).all():
        return False
    else:
        return True
    
if plotting_available:
    def plot_cosmonuc(results, logscale = True, *args, **kwargs):
        """Plots 2D histogram of cosmogenic nuclides. Pass the planetoparse
        instance to plot."""
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
            c = log10(results.cosmonuc.data[:, 4])
        else:
            c = results.cosmonuc.data[:, 4]
        plt.scatter(results.cosmonuc.data[:, 0] + .5, 
                    results.cosmonuc.data[:, 2] + .5, c = c, s = s, 
                    marker = marker, lw = lw, *args, **kwargs)
        plt.xlabel(results.cosmonuc.params['Xaxis'])
        plt.ylabel(results.cosmonuc.params['Yaxis'])
        plt.xlim([0, 25])
        plt.xticks(arange(0, 26, 2))
        plt.ylim([0, 25])
        plt.yticks(arange(0, 26, 2))
        fig = plt.gcf()
        ax = fig.get_axes()
        title, units = __parse_title(results.cosmonuc.params['Title'])
        #remove colorbar, if already present
        if len(ax) == 2:
            if ax[0].get_title() == title and ax[1].get_title() == '':
                fig.delaxes(ax[1])
                fig.subplots_adjust()
        cbar = plt.colorbar()
        plt.title(title)
        if logscale:
            units = 'log ' + units
        cbar.set_label(units)
        plt.show(block = False)
        return
    
if plotting_available:
    def plot_detector_levels(fluxhists, plot_only = [], dont_plot = []):
        """Plots the detector levels of the given flux histograms into the 
        current plot."""
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
            altitude, alt_unit = \
                fluxhists[detector].params['Altitude'].split(' ')
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
    
def combine_histograms(hists, scale_by = 1.):
    """Combines the different histograms for multiple planetoparse instances.
    WARNING: Only the primaries count, cosmonuc 2D histograms, primary 
    particle fluxes and up/down flux histograms are combined!"""
    if len(hists) > 1:
        res = planetoparse()
        for addthis in hists:
            #check and set normalisation
            if not res.normalisation == '':
                if not res.normalisation == addthis.normalisation:
                    print 'WARNING: Skipping result, different normalisation\
                        methods'
                    print '\t' + addthis.normalisation + ' vs. ' +\
                        res.normalisation + ' as expected'
                    continue
            else:
                res.normalisation += addthis.normalisation
            #combine primary counts
            res.primaries += addthis.primaries
            #combine number of events counts
            if 'nb_of_events' in res.params:
                res.params['nb_of_events'] += addthis.params['nb_of_events']
            else:
                res.params['nb_of_events'] = float64(0.) +\
                    addthis.params['nb_of_events']
            #combine cosmogenic nuclides
            if not res.cosmonuc is None:
                res.cosmonuc.data[:, 4] += addthis.cosmonuc.data[:, 4]
            else:
                res.cosmonuc = histdata(copyhist = addthis.cosmonuc)
            #combine upward fluxes
            for particle in addthis.flux_up:
                for detector in addthis.flux_up[particle]:
                    if not particle in res.flux_up:
                        res.flux_up[particle] = {}
                    if detector in res.flux_up[particle]:
                        res.flux_up[particle][detector] = \
                            __combine_single_hists(
                                res.flux_up[particle][detector], 
                                addthis.flux_up[particle][detector], 
                                scale_by = scale_by)
                    else:
                        res.flux_up[particle][detector] = \
                            histdata(copyhist = \
                                addthis.flux_up[particle][detector])
            #combine downward fluxes
            for particle in addthis.flux_down:
                for detector in addthis.flux_down[particle]:
                    if not particle in res.flux_down:
                        res.flux_down[particle] = {}
                    if detector in res.flux_down[particle]:
                        res.flux_down[particle][detector] = \
                            __combine_single_hists(
                                res.flux_down[particle][detector], 
                                addthis.flux_down[particle][detector], 
                                scale_by = scale_by)
                    else:
                        res.flux_down[particle][detector] = \
                            histdata(copyhist = \
                                addthis.flux_down[particle][detector])
            #combine primary fluxes
            for particle in addthis.primhists:
                if not particle in res.primhists:
                    res.primhists[particle] = []
                    for hist in addthis.primhists[particle]:
                        res.primhists[particle].append(
                            histdata(copyhist = hist))
                else:
                    for i in arange(0, len(addthis.primhists[particle])):
                        added = False
                        for j in arange(0, len(res.primhists[particle])):
                            addtitle = __parse_title(
                                addthis.primhists\
                                    [particle][i].params['Title'])[0]
                            restitle = __parse_title(res.primhists\
                                [particle][j].params['Title'])[0]
                            if addtitle == restitle:
                                res.edep_atmo[j] = __combine_single_hists(
                                    res.primhists[particle][j],
                                    addthis.primhists[particle][i], 
                                    scale_by = scale_by)
                                added = True
                        if not added:
                            res.edep_atmo.append(histdata(
                                copyhist = addthis.primhists[particle][i]))
                            added = True
            #combine 1d hist list
            for i in arange(0, len(addthis.hists1d)):
                added = False
                for j in arange(0, len(res.hists1d)):
                    addtitle = __parse_title(
                        addthis.hists1d[i].params['Title'])[0]
                    restitle = __parse_title(
                        res.hists1d[j].params['Title'])[0]
                    if addtitle == restitle and addthis.hists1d[i].detector\
                        == res.hists1d[j].detector:
                        res.hists2d[j] = __combine_single_hists(
                            res.hists1d[j], addthis.hists1d[i], 
                            scale_by = scale_by)
                        added = True
                if not added:
                    res.hists2d.append(histdata(
                        copyhist = addthis.hists2d[i]))
                    added = True
            #combine 2d hist list
            for i in arange(0, len(addthis.hists2d)):
                added = False
                for j in arange(0, len(res.hists2d)):
                    addtitle = __parse_title(
                        addthis.hists2d[i].params['Title'])[0]
                    restitle = __parse_title(
                        res.hists2d[j].params['Title'])[0]
                    if addtitle == restitle and addthis.hists2d[i].detector\
                        == res.hists2d[j].detector:
                        res.hists2d[j] = __combine_single_hists(
                            res.hists2d[j], addthis.hists2d[i], 
                            scale_by = scale_by)
                        added = True
                if not added:
                    res.hists2d.append(histdata(
                        copyhist = addthis.hists2d[i]))
                    added = True
            #combine 2d upward fluxes
            for particle in addthis.flux2d_up:
                for detector in addthis.flux2d_up[particle]:
                    if not particle in res.flux2d_up:
                        res.flux2d_up[particle] = {}
                    if detector in res.flux2d_up[particle]:
                        for i in arange(0, len(
                            addthis.flux2d_up[particle][detector])):
                            added = False
                            for j in arange(0, len(res.flux2d_up\
                                [particle][detector])):
                                addtitle = __parse_title(addthis.flux2d_up\
                                    [particle][detector][i].params['Title'])\
                                        [0]
                                restitle = __parse_title(
                                    res.flux2d_up[particle][detector][j].\
                                        params['Title'])[0]
                                if addtitle == restitle:
                                    res.flux2d_up[particle][detector][j] = \
                                        __combine_single_hists(
                                            res.flux2d_up[particle]\
                                                [detector][j], 
                                        addthis.flux2d_up[particle]\
                                            [detector][i], scale_by = scale_by)
                                    added = True
                            if not added:
                                res.flux2d_up[particle][detector].append(
                                    histdata(copyhist = \
                                        addthis.flux2d_up[particle]\
                                        [detector][i]))
                                added = True
                    else:
                        res.flux2d_up[particle][detector] = []
                        for hist in addthis.flux2d_up[particle][detector]:
                            res.flux2d_up[particle][detector].append(
                                histdata(copyhist = hist))
            #combine 2d downward fluxes
            for particle in addthis.flux2d_down:
                for detector in addthis.flux2d_down[particle]:
                    if not particle in res.flux2d_down:
                        res.flux2d_down[particle] = {}
                    if detector in res.flux2d_down[particle]:
                        for i in arange(0, len(addthis.flux2d_down[particle]\
                            [detector])):
                            added = False
                            for j in arange(0, len(res.flux2d_down[particle]\
                                [detector])):
                                addtitle = __parse_title(addthis.flux2d_down\
                                    [particle][detector][i].params['Title'])\
                                        [0]
                                restitle = __parse_title(
                                    res.flux2d_down[particle][detector][j].\
                                        params['Title'])[0]
                                if addtitle == restitle:
                                    res.flux2d_down[particle][detector][j]\
																			= __combine_single_hists(
                                            res.flux2d_down[particle]\
                                                [detector][j], addthis.\
                                                flux2d_down[particle]\
                                                [detector][i], 
                                                scale_by = scale_by)
                                    added = True
                            if not added:
                                res.flux2d_down[particle][detector].append(
                                    histdata(copyhist = addthis.flux2d_down\
                                        [particle][detector][i]))
                                added = True
                    else:
                        res.flux2d_down[particle][detector] = []
                        for hist in addthis.flux2d_down[particle][detector]:
                            res.flux2d_down[particle][detector].append(
                                histdata(copyhist = hist))
            #combine atmospheric edep hists
            for i in arange(0, len(addthis.edep_atmo)):
                added = False
                for j in arange(0, len(res.edep_atmo)):
                    addtitle = __parse_title(
                        addthis.edep_atmo[i].params['Title'])[0]
                    restitle = __parse_title(
                        res.edep_atmo[j].params['Title'])[0]
                    if addtitle == restitle:
                        res.edep_atmo[j] = __combine_single_hists(
                            res.edep_atmo[j], addthis.edep_atmo[i], 
                            scale_by = scale_by)
                        added = True
                if not added:
                    res.edep_atmo.append(
                        histdata(copyhist = addthis.edep_atmo[i]))
                    added = True
            #combine soil edep hists
            for i in arange(0, len(addthis.edep_soil)):
                added = False
                for j in arange(0, len(res.edep_soil)):
                    addtitle = __parse_title(
                        addthis.edep_soil[i].params['Title'])[0]
                    restitle = __parse_title(
                        res.edep_soil[j].params['Title'])[0]
                    if addtitle == restitle:
                        res.edep_soil[j] = __combine_single_hists(
                            res.edep_soil[j], addthis.edep_soil[i], 
                            scale_by = scale_by)
                        added = True
                if not added:
                    res.edep_soil.append(histdata(
                        copyhist = addthis.edep_soil[i]))
                    added = True
            #combine atmo ionization hists
            for i in arange(0, len(addthis.edep_ionization)):
                added = False
                for j in arange(0, len(res.edep_ionization)):
                    addtitle = __parse_title(
                        addthis.edep_ionization[i].params['Title'])[0]
                    restitle = __parse_title(
                        res.edep_ionization[j].params['Title'])[0]
                    if addtitle == restitle:
                        res.edep_ionization[j] = __combine_single_hists(
                            res.edep_ionization[j], addthis.edep_ionization[i], 
                            scale_by = scale_by)
                        added = True
                if not added:
                    res.edep_ionization.append(histdata(
                        copyhist = addthis.edep_ionization[i]))
                    added = True
        return res
    else:
        return hists[0]
        
def __combine_single_hists(hist1, hist2, scale_by = 1.):
    """Combine two histograms into one."""
    from planetoparse import histdata
    if not hist1.type == hist2.type:
        print 'ERROR: Unable to combine histograms, incompatible dimensions'
        return hist1
    elif hist1.type == hist2.type == 'Histogram1D':
        res = histdata(copyhist = hist1)
        if not (res.data[:, 0] == hist2.data[:, 0]).all() or not \
            (res.data[:, 1] == hist2.data[:, 1]).all() or not \
            (res.data[:, 2] == hist2.data[:, 2]).all():
            print 'WARNING: Histogram binning is different, interpolating.'
            xunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
            yunits1 = __normalize_units(__parse_title(res.params['Title'])[1])
            xunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
            yunits2 = __normalize_units(__parse_title(hist2.params['Title'])[1])
            if not (xunits1 == xunits2 and yunits1 == yunits2):
                print 'ERROR: Unable to combine histograms, units mismatch.'
                return hist1
            from scipy.interpolate import interp1d
            interpolator = interp1d(hist2.data[:, 2], hist2.data[:, 3],
                                    bounds_error = False, fill_value = 0.)
            error_interpolator = interp1d(hist2.data[:, 2], hist2.data[:, 4],
                                          bounds_error = False, fill_value = 0.)
            res.data[:, 3] += interpolator(res.data[:, 2])
            res.data[:, 3] /= scale_by
            res.data[:, 4] = sqrt(res.data [:, 4]**2 + 
                                  error_interpolator(res.data[:, 2])**2)/scale_by
            return res            
        xunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
        yunits1 = __normalize_units(__parse_title(res.params['Title'])[1])
        xunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
        yunits2 = __normalize_units(__parse_title(hist2.params['Title'])[1])
        if not (xunits1 == xunits2 and yunits1 == yunits2):
            print 'ERROR: Unable to combine histograms, units mismatch.'
            return hist1
        res.data[:, 3] += hist2.data[:, 3]
        res.data[:, 3] /= scale_by
        res.data[:, 4] = sqrt(res.data [:, 4]**2 + hist2.data[:, 4]**2)/scale_by
        return res
    elif hist1.type == hist2.type == 'Histogram2D':
        res = histdata(copyhist = hist1)
        if not (res.data[:, 0] == hist2.data[:, 0]).all() or not \
            (res.data[:, 1] == hist2.data[:, 1]).all() or not \
            (res.data[:, 2] == hist2.data[:, 2]).all() or not \
            (res.data[:, 3] == hist2.data[:, 3]).all():
            print 'WARNING: Histogram binning is different, interpolating.'
            xunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
            yunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
            zunits1 = __normalize_units(__parse_title(res.params['Title'])[1])
            xunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
            yunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
            zunits2 = __normalize_units(__parse_title(hist2.params['Title'])[1])
            if not (xunits1 == xunits2 and yunits1 == yunits2 and\
                zunits1 == zunits2):
                print 'ERROR: Unable to combine histograms, units mismatch.'
                return hist1
            from scipy.interpolate import interp1d
            interpolator = interp2d(hist2.data[:, 1], hist2.data[:, 3], 
                                    hist2.data[:, 4], bounds_error = False, 
                                    fill_value = 0.)
            error_interpolator = interp2d(hist2.data[:, 1], hist2.data[:, 3], 
                                          hist2.data[:, 5], bounds_error = False, 
                                          fill_value = 0.)
            res.data[:, 4] += interpolator(res.data[:, 1], res.data[:, 3])
            res.data[:, 4] /= scale_by
            res.data[:, 5] = sqrt(res.data [:, 5]**2 + 
                                  error_interpolator(res.data[:, 1], 
                                               res.data[:, 3])**2)/scale_by
            return res            
        xunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
        yunits1 = __normalize_units(__parse_title(res.params['Xaxis'])[1])
        zunits1 = __normalize_units(__parse_title(res.params['Title'])[1])
        xunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
        yunits2 = __normalize_units(__parse_title(hist2.params['Xaxis'])[1])
        zunits2 = __normalize_units(__parse_title(hist2.params['Title'])[1])
        if not (xunits1 == xunits2 and yunits1 == yunits2 and\
            zunits1 == zunits2):
            print 'ERROR: Unable to combine histograms, units mismatch.'
            return hist1
        res.data[:, 4] += hist2.data[:, 4]
        res.data[:, 4] /= scale_by
        res.data[:, 5] = sqrt(res.data [:, 5]**2 + hist2.data[:, 5]**2)/scale_by
        return res
        
def add_histograms(hists, scale_by = 1.):
    """Add single histograms."""
    if len(hists) > 1:
        res = histdata(copyhist = hists[0])
        for hist in hists[1:]:
            res = __combine_single_hists(res, hist, scale_by = scale_by)
        return res
    else:
        return hists[0]
            
def get_normalization_factors(results):
    """Get a dictionary with the normalization factors for all particles in 
    the flux_up and flux_down sections of the given results."""
    res = {}
    for flux_list in [results.flux_down, results.flux_up]:
        for particle in flux_list:
            #this is a completely arbitrary choice
            detector = flux_list[particle].keys()[-1] 
            if 'normalisation_factor' in flux_list\
                [particle][detector].params:
                factor = flux_list[particle][detector].params\
                    ['normalisation_factor']
            else:
                factor = nan
            if not particle in res:
                res[particle] = factor
    return res
    
def print_list_info(histlist, indices = None):
    """Prints information on histograms in the given histogram list. 
    Optional second argument selects the indices to print."""
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
        message += '{index: >3d}'.format(index = index) + \
            ': {particle: <10s}'.format(particle = histlist[index].particle)\
            + ', detector {detector: >3d}'.format(
                detector = histlist[index].detector) + ': '
        message += histlist[index].params['Title'] + '\n'
    print message[:-1]
    return

def integrate_fluxhist(histogram, limits = None):
    """Integrates the flux stored in the given 1D histogram. Pass a tuple 
    (Emin, Emax) as second argument to set the integration range."""
    if isinstance(histogram, histdata):
        if not limits is None:
            lims_min = amin(limits)
            lims_max = amax(limits)
            mask = (histogram.data[:, 0] >= lims_min) * \
                (histogram.data[:, 1] <= lims_max)
        else:
            mask = ones(len(histogram.data), dtype = bool)
        return sum(histogram.data[:, 3][mask])
    elif type(histogram) == dict:
        if not limits is None:
            lims_min = amin(limits)
            lims_max = amax(limits)
            mask = (histogram['x'] >= lims_min) * \
                (histogram['x'] <= lims_max)
        else:
            mask = ones(len(histogram['x']), dtype = bool)
        return sum(histogram['y'][mask])
        
def integrate_2d_fluxhist(histogram, xlimits = None, ylimits = None):
    """Integrates the flux stored in the given 2D histogram. Pass a tuple 
    (Emin, Emax) as xlimits or ylimits keyword arguments to set the 
    integration range."""
    if not xlimits is None:
        xlims_min = amin(xlimits)
        xlims_max = amax(xlimits)
        xmask = (histogram.data[:, 0] >= xlims_min) * \
            (histogram.data[:, 1] <= xlims_max)
    else:
        xmask = ones(len(histogram.data), dtype = bool)
    if not ylimits is None:
        ylims_min = amin(ylimits)
        ylims_max = amax(ylimits)
        ymask = (histogram.data[:, 2] >= ylims_min) * \
            (histogram.data[:, 3] <= ylims_max)
    else:
        ymask = ones(len(histogram.data), dtype = bool)
    mask = xmask * ymask
    return sum(histogram.data[:, 4][mask])
    
def project_data(histogram, axis = 'x', xlimits = None, ylimits = None):
    """Project data in a 2D histogram onto either the X or the Y axis 
    and return the resulting 1D histogram."""
    if not axis in ['x', 'y']:
        print 'ERROR: Invalid axis selection'
        return None
    if not histogram.type == 'Histogram2D':
        print 'ERROR: Cannot project 1D histograms'
        return
    #apply selected limits
    if not xlimits is None:
        xlims_min = amin(xlimits)
        xlims_max = amax(xlimits)
        xmask = (histogram.data[:, 0] >= xlims_min) * \
            (histogram.data[:, 1] <= xlims_max)
    else:
        xmask = ones(len(histogram.data), dtype = bool)
    if not ylimits is None:
        ylims_min = amin(ylimits)
        ylims_max = amax(ylimits)
        ymask = (histogram.data[:, 2] >= ylims_min) * \
            (histogram.data[:, 3] <= ylims_max)
    else:
        ymask = ones(len(histogram.data), dtype = bool)
    mask = xmask * ymask
    data = empty_like(histogram.data)
    data[:] = histogram.data
    data = data[mask]
    #find data binning for given axis
    if axis == 'x':
        index_lower = 0
        index_upper = 1
    else:
        index_lower = 2
        index_upper = 3
    binslower = []
    binsupper = []
    for line in data:
        if not line[index_lower] in binslower:
            binslower.append(line[index_lower])
        if not line[index_upper] in binsupper:
            binsupper.append(line[index_upper])
    binslower = array(sort(binslower))
    binsupper = array(sort(binsupper))
    #build data array
    bindata = zeros(len(binslower))
    binerrors = zeros(len(binslower))
    for line in data:
        index = -1
        for i in arange(0, len(bindata)):
            if line[index_lower] == binslower[i] and \
                line[index_upper] == binsupper[i]:
                index = i
                break
        if not index == -1:
            bindata[index] += line[4]
            binerrors[index] = sqrt(binerrors[index]**2 + line[5]**2)
        else:
            print 'WARNING: Skipping data point, bin not found'
    binmiddles = binslower + (binsupper - binslower) / 2
    resdata = column_stack((binslower, binsupper, 
                            binmiddles, bindata, binerrors))
    #build result histdata object
    res = histdata()
    res.data = resdata
    #build and set parameter dict
    for entry in histogram.params:
        if 'axis' in entry:
            if entry == axis.upper() + 'axis':
                res.params['Xaxis'] = '' + histogram.params[entry]
        else:
            if type(histogram.params[entry]) == str:
                res.params[entry] = '' + histogram.params[entry]
            else:
                res.params[entry] = float64(0.) + histogram.params[entry]
    #build and set title string
    title = ''
    #add information about limits applied to data
    limits_string = ''
    if not xlimits is None:
        limits_string += ', '
        tmp1, tmp2 = __parse_title(histogram.params['Xaxis'])
        if tmp1 is None:
            xlabel = tmp2
        else:
            xlabel = tmp1
        limits_string += str(amin(xlimits)) + '<=' + xlabel + '<=' + \
            str(amax(xlimits))
    if not ylimits is None:
        limits_string += ', '
        tmp1, tmp2 = __parse_title(histogram.params['Yaxis'])
        if tmp1 is None:
            ylabel = tmp2
        else:
            ylabel = tmp1
        limits_string += str(amin(ylimits)) + '<=' + ylabel + '<=' + \
            str(amax(ylimits))
    tmp1, tmp2 = __parse_title(histogram.params[axis.upper() + 'axis'])
    if tmp1 is None:
        axislabel = tmp2
    else:
        axislabel = tmp1    
    title += 'Flux vs ' + axislabel + ' of ' + histogram.particle 
    title += limits_string
    title += '[' + __parse_title(histogram.params['Title'])[1] + ']'
    res.params['Title'] = title
    #set remaining parameters
    res.detector = 0 + histogram.detector
    res.particle = '' + histogram.particle
    res.type = 'Histogram1D'
    return res

def enumerate_hist_list(list):
    """Print list index, particle type and detector for a given list of 
    histograms."""
    for i in xrange(0, len(list)):
        print '{i:3d}: det. {det:2d} {part:s}'.format(i=i, 
                                                      part=list[i].particle, 
                                                      det=list[i].detector)
                                                      
def hist2dose(hist, Z_target=3.33333, A_target=6., rho_target=1e3,
              I=75., relativistic=True, **kwargs):
    """Compute dose rate from an energy spectrum. Defaults to dose in H2O."""
    from bethebloch import bethebloch
    from scipy.constants import eV, m_u, c
    from numpy import isnan
    doserates = []
    doserates_delta = []
    for (E_low, E_high, E, flux, flux_delta) in hist.data:
        E = E*1e6*eV #E in Joule
        E_low = E_low*1e6*eV
        E_high = E_high*1e6*eV
        E_delta = (E_high - E_low)/2
        flux = flux*1e4 #flux in 1/m^2/s
        flux_delta = flux_delta*1e4
        (Z_proj, A_proj) = __get_Z_A(hist.particle)
        m = A_proj*m_u #m in kg
        if not relativistic:
            v_proj = sqrt(2*E/m) #V in m/s
            v_proj_delta = sqrt(1/2/m/E_delta)
        else:
            if E < m*c**2:
                E = 0.
            v_proj = c*sqrt(1-m**2*c**4/(E - m*c**2)**2)
            v_proj_delta = 0.
        if isnan(v_proj):
            print 'v is nan'
            print '\tE =', E
        if v_proj > c:
            print 'v > c, constraining'
            print '\tE =', E
            v_proj = .99999*c
        dEdx = bethebloch(v_proj, Z_proj, Z_target, A_target, rho_target,
                          I=I, **kwargs) #dE/dx in J/m
        doserates.append(flux*dEdx/rho_target)
        #FIXME: At some point, we should compute the delta of the dose rate...
    return sum(doserates)

def __get_Z_A(part_name):
    """Get particle Z and A from the Geant4 particle name"""
    from re import match
    part_names = {'Ac227': (89, 227), 'Ag107': (47, 107), 'Al27': (13, 27),
                  'Am243': (95, 243), 'Ar40': (18, 40), 'As75': (33, 75),
                  'At210': (85, 210), 'Au197': (79, 197), 'B11': (5, 11),
                  'Ba138': (56, 138), 'Be9': (4, 9), 'Bi209': (83, 209),
                  'Br79': (35, 79), 'C12': (6, 12), 'Ca40': (20, 40),
                  'Cd114': (48, 114), 'Ce140': (58, 140), 'Cl35': (17, 35),
                  'Cm247': (96, 247), 'Co59': (27, 59), 'Cr52': (24, 52),
                  'Cs133': (55, 133), 'Cu63': (29, 63), 'Dy164': (66, 164),
                  'Er166': (68, 166), 'Eu153': (63, 153), 'F19': (9, 19),
                  'Fe56': (26, 56), 'Fr223': (87, 223), 'Ga69': (31, 69),
                  'Gd158': (64, 158), 'Ge74': (32, 74), 'Hf180': (72, 180),
                  'Hg202': (80, 202), 'Ho165': (67, 165), 'I127': (53, 127),
                  'In115': (49, 115), 'Ir193': (77, 193), 'K39': (19, 39),
                  'Kr84': (36, 84), 'La139': (57, 139), 'Li7': (3, 7),
                  'Lu175': (71, 175), 'Mg24': (12, 24), 'Mn55': (25, 55),
                  'Mo98': (42, 98), 'N14': (7, 14), 'Na23': (11, 23),
                  'Nb93': (41, 93), 'Nd144': (60, 144), 'Ne20': (10, 20),
                  'Ni58': (28, 58), 'Np237': (93, 237), 'O16': (8, 16),
                  'Os192': (76, 192), 'P31': (15, 31), 'Pa231': (91, 231),
                  'Pb208': (82, 208), 'Pd106': (46, 106), 'Pm145': (61, 145),
                  'Po209': (84, 209), 'Pr141': (59, 141), 'Pt195': (78, 195),
                  'Pu244': (94, 244), 'Ra226': (88, 226), 'Rb58': (37, 58),
                  'Re187': (75, 187), 'Rh103': (45, 103), 'Rn222': (86, 222),
                  'Ru102': (44, 102), 'S32': (16, 32), 'Sb121': (51, 121),
                  'Sc45': (21, 45), 'Se80': (34, 80), 'Si28': (14, 28), 
                  'Sm152': (62, 152), 'Sn120': (50, 120), 'Sr88': (38, 88),
                  'Ta181': (73, 181), 'Tb159': (65, 159), 'Tc98': (43, 98),
                  'Te130': (52, 130), 'Th232': (90, 232), 'Ti48': (22, 48),
                  'Tl205': (81, 205), 'Tm169': (69, 169), 'U238': (92, 238),
                  'V51': (23, 51), 'W184': (74, 184), 'Xe132': (54, 132),
                  'Y89': (39, 89), 'Yb174': (70, 174), 'Zn64': (30, 64),
                  'Zr90': (40, 90), 'alpha': (2, 4), 'proton': (1, 1),
                  'e+': (1, 5.485e-4), 'e-': (-1, 5.485e-4), 
                  'mu+': (1, 1.134e-1), 'mu-': (-1, 1.134e-1), 
                  'pi+': (1, 1.498e-1), 'pi-': (-1, 1.498e-1)}
    part_name = part_name.replace('[0.0]', '') #Fix stupid Geant4.9 naming
    if not part_name in part_names:
        #Try to find the element from combined high-Z isotope hists:
        for name in part_names:
            res = match('([A-Za-z]+)[0-9]+', name)
            if not res is None:
                if res.group(1) == part_name:
                    return part_names[name]
        #everything failed, return nan
        return (nan, nan)
    else:
        return part_names[part_name]
