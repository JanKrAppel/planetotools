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
    
def plot_data_overview(mcddata, title = False, legend = False, *args, **kwargs):
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
    if legend:
        make_overview_legend()
    return
    
def make_overview_legend():
    """Plot a legend into all four data overview subplots."""
    plt.subplot(221)
    plt.legend(loc = 'best')
    plt.subplot(222)
    plt.legend(loc = 'best')
    plt.subplot(223)
    plt.legend(loc = 'best')
    plt.subplot(224)
    plt.legend(loc = 'best')
    return
    
def __get_daily_mean_dustdepth(date, dustscen, lon, lat):
    day, month, year, hour, minute, second = date
    hours = arange(0, 23, 1)
    dust = []
    ls = []
    loct = []
    for hour in hours:
        (ierr, xdate) = julian(month, day, year, hour, 0, 0)
        (pres, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr) = call_mcd(3, 0, lon, lat, 1, 0, xdate, 0, MCD_DATA_DIR, dustscen, 4, 10, 0, ones(100))
        if extvar[35] == -999.0:
            res = nan
        else:
            res = extvar[35]
        dust.append(res)
        ls.append(extvar[4])
        loct.append(extvar[5])
    return mean(dust), amax(dust), amin(dust), mean(ls), mean(loct)

def plot_yearly_dust(year, dust, lon = 0., lat = 0., plot_only = None, show_marstime = True, *args, **kwargs):
    """Plots the daily mean, min and max values for the dust optical depth for the given year and dust scenario."""
    #parse dust scenario, if needed
    if type(dust) == str:
        if dust in dust_scenarios:
            dust = dust_scenarios[dust]
        else:
            print 'ERROR: Invalid dust scenario selected.'
            return
    dustcurve = []
    maxcurve = []
    mincurve = []
    daymonth = []
    marstime = []
    for month in arange(1, 13, 1):
        for day in arange(1, 31, 1):
            dustmean, dustmax, dustmin, ls, loct = __get_daily_mean_dustdepth((day, month, year, 12, 00, 00), dust, lon = lon, lat = lat)
            daymonth.append(str(day) + '/' + str(month))
            marstime.append('Ls ' + str(around(ls, 1)))
            dustcurve.append(dustmean)
            maxcurve.append(dustmax)
            mincurve.append(dustmin)
    if plot_only is None:
        if 'label' in kwargs:
            label = kwargs['label']
            kwargs.pop('label')
        else:
            label = 'mean'
        plt.plot(dustcurve, 'bo', label = label, *args, **kwargs)
        if 'label' in kwargs:
            label = kwargs['label']
            kwargs.pop('label')
        else:
            label = 'max'
        plt.plot(maxcurve, 'ro', label = label, *args, **kwargs)
        if 'label' in kwargs:
            label = kwargs['label']
            kwargs.pop('label')
        else:
            label = 'min'
        plt.plot(mincurve, 'go', label = label, *args, **kwargs)
    else:
        if 'mean' in plot_only:
            if 'label' in kwargs:
                label = kwargs['label']
                kwargs.pop('label')
            else:
                label = 'mean'
            plt.plot(dustcurve, 'o', label = label, *args, **kwargs)
        if 'max' in plot_only:
            if 'label' in kwargs:
                label = kwargs['label']
                kwargs.pop('label')
            else:
                label = 'max'
            plt.plot(maxcurve, 'o', label = label, *args, **kwargs)
        if 'min' in plot_only:
            if 'label' in kwargs:
                label = kwargs['label']
                kwargs.pop('label')
            else:
                label = 'min'
            plt.plot(mincurve, 'o', label = label, *args, **kwargs)
    plt.legend(loc = 'best')
    #Earth time xaxis:
    locs, labels = plt.xticks()
    daymonth = array(daymonth)
    locs = locs[locs<= len(daymonth)]
    newlabels = daymonth[int64(locs)]
    plt.xticks(locs, newlabels)
    plt.xlabel('Day/Month')
    #Mars time xaxis:
    if show_marstime:
        marstime = array(marstime)
        newlabels = marstime[int64(locs)]
        plt.gca().twiny() #get shared y axis -> second x axis
        plt.xticks(locs, newlabels)
        plt.xlabel('Mars solar longitude')
    #set ylabel and title
    plt.ylabel('Dust optical depth')
    for value in dust_scenarios:
        if dust_scenarios[value] == dust:
            dustscen = value
    if show_marstime:
        plt.suptitle('Dust optical depth variations for ' + str(year) + ', scenario ' + dustscen)
    else:
        plt.title('Dust optical depth variations for ' + str(year) + ', scenario ' + dustscen)
    return    

