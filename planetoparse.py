#!/usr/bin/env python
"""Contains two class definitions:
histdata contains histogram data and metadata as well as some functionality
for scaling, saving and loading histograms.
planetoparse is the central planetotools data class, it parses Planetocosmics
ASCII output files, stores the histograms in histdata objects in a structured
way and provides scaling, saving and loading functionality."""

from numpy import *
import re
import cPickle

DELIMITER = '//////////////////////////////////////////////////\n'
PARAMS_PATTERN = '(\w*)\s*:\s*(.*)\n'
DEFAULT_SORT_CONFIG = {'flux_up': ['UpHist', '2'], 
                       'flux_down': ['DownHist', '1'], 
                       'flux2d_up': ['UpPosHist', '2'],
                       'flux2d_down': ['DownPosHist', '1'], }

####################
#histdata class definition
####################

class histdata:
    """Provides histogram data to plotting routines. Allows scaling and 
    unscaling per nucleus and per sterad, as well as saving and loading 
    to and from Planetocosmics flux data definition files."""
    
    def __init__(self, copyhist = None):
        self.type = ''
        self.title = ''
        self.params = {}
        self.data = ''
        self.particle = ''
        self.detector = -1
        self.scaled_per_nuc = False
        self.nuc_weight = 0
        self.scaled_per_sterad = False
        if not copyhist is None:
            self.type = '' + copyhist.type
            self.title = '' + copyhist.title
            for param in copyhist.params:
                if isinstance(copyhist.params[param], float64):    
                    self.params[param] = float64(0.) + copyhist.params[param]
                else:
                    self.params[param] = '' + copyhist.params[param]
            self.data = zeros(copyhist.data.shape)
            for i in arange(0, len(copyhist.data), 1):
                self.data[i] += copyhist.data[i]
            self.particle = '' + copyhist.particle
            self.detector = 0 + copyhist.detector
            self.nuc_weight = 0 + copyhist.nuc_weight
            if copyhist.scaled_per_nuc:
                self.scaled_per_nuc = True
            if copyhist.scaled_per_sterad:
                self.scaled_per_sterad = True

    def scale_per_nuc(self, weight = None):
        """Scale histogram to Energy/nuc. Pass the nucleus weight in amu."""
        if weight is None and self.nuc_weight == 0:
            return False
        if not self.scaled_per_nuc:
            if not weight is None:
                self.nuc_weight = weight
            self.data[:, :3] /= self.nuc_weight
            self.data[:, 3:] *= self.nuc_weight
            titleparse = re.match('(.*)\s*\[(.*)\]', self.params['Xaxis'])
            self.params['Xaxis'] = \
                titleparse.group(1) + '[' +  titleparse.group(2) + '/nuc]'
            self.scaled_per_nuc = True
            return True
        else:
            return False
            
    def unscale_per_nuc(self):
        """Remove energy scaling per nucleon from histogram."""
        if self.scaled_per_nuc:
            self.data[:, :3] *= self.nuc_weight
            self.data[:, 3:] /= self.nuc_weight
            self.params['Xaxis'] = re.sub('/nuc', '', self.params['Xaxis'])
            self.scaled_per_nuc = False
            return True
        else:
            return False
                        
    def scale_per_sterad(self):
        """Scale histogram to Energy/sr."""
        if not self.scaled_per_sterad:
            self.data[:, 3] /= 2*pi
            self.data[:, 4] /= 2*pi
            titleparse = re.match('(.*)\s*\[(.*)\]', self.params['Xaxis'])
            self.params['Xaxis'] = \
                titleparse.group(1) + '[' +  titleparse.group(2) + '/sr]'
            self.scaled_per_sterad = True
            return True
        else:
            return False
            
    def unscale_per_sterad(self):
        """Remove energy scaling per sr from histogram."""
        if self.scaled_per_sterad:
            self.data[:, 3] *= 2*pi
            self.data[:, 4] *= 2*pi
            self.params['Xaxis'] = re.sub('/sr', '', self.params['Xaxis'])
            self.scaled_per_sterad = False
            return True
        else:
            return False
            
    def isempty(self):
        """Returns true if the histogram is all-zero, false if it is not.
        None is returned when no valid histogram is loaded."""
        if self.type == 'Histogram1D':
            return self.__histogram_1d_empty()
        elif self.type == 'Histogram2D':
            return self.__histogram_2d_empty()
        else:
            return None
            
    def __histogram_2d_empty(self):
        """Checks if a 2D histogram is empty, returns boolean"""
        return (self.data[:, 4] == zeros(len(self.data[:, 4]))).all()

    def __histogram_1d_empty(self):
        """Checks if a 1D histogram is empty, returns boolean"""
        return (self.data[:, 3] == zeros(len(self.data[:, 3]))).all()

    def save_as_flux(self, filename):
        """Save the histogram in a format that can be read in Planetocosmics
        for primary flux definition. Pass the output filename."""
        if self.type == 'Histogram1D':
            def get_unit(title):
                """Parses title information and returns the units."""
                titleparse = re.match('(.*)\s*\[(.*)\]', title)
                if not titleparse is None:
                    return titleparse.group(2)
                else:
                    return ''
            eunit = get_unit(self.params['Xaxis'])
            fluxunit = get_unit(self.params['Title'])
            fluxunit = re.sub('nb particles', '#', fluxunit)
            particle = re.sub('primary ', '', self.particle)
            outfile = open(filename, 'w')
            outfile.write('\\definition\n')
            outfile.write('\\energy_unit{' + eunit + '}\n')
            outfile.write('\\flux_unit{' + fluxunit + '}\n')
            outfile.write('\\particle{' + particle + '}\n')
            outfile.write('\\interpolation{log}\n')
            outfile.write('\\data\n')
            for i in arange(0, len(self.data), 1):
                outfile.write(str(self.data[i, 2]) + '\t' \
                    + str(self.data[i, 3]) + '\n')
            outfile.close()
            return
        else:
            print 'ERROR: Can only save 1D histogram data as flux \
                definition.'
            return
            
    def save_data(self, filename):
        """Save only the data array to disk."""
        savetxt(filename, self.data)
        return
            
    def __parse_params(self, line):
        """Parses parameter information. Returns name, value pair if found,
        None, None pair if not found."""
        parsed = re.match('\\\\(\w.*)\{(.*)\}', line)
        if parsed:
            name = parsed.group(1)
            try:
                value = float64(parsed.group(2))
            except ValueError:
                value = parsed.group(2)
            return name, value
        else:
            return None, None

    def load_from_flux(self, filename):
        """Load the histogram data from a file containing Planetocosmics
        primary flux definitions."""
        infile = open(filename, 'r')
        line = infile.readline()
        tmpparams = {}
        tmpdata = []
        while not line == '':
            if '\\definition' in line:
                line = infile.readline()
                while not '\\data' in line:
                    name, value = self.__parse_params(line)
                    if not name is None:
                        tmpparams[name] = value
                    line = infile.readline()
            elif '\\data' in line:
                line = infile.readline()
                while not line[:-1] == '':
                    tmpdata.append(array(line.split(), dtype = float64))
                    line = infile.readline()
                tmpdata = array(tmpdata)
        infile.close()
        #set type, particle and title
        self.type = 'Histogram1D'
        self.title = '/PRIMARY/' + tmpparams['particle'] + '/1'
        self.particle = 'primary ' + tmpparams['particle']
        #set other params
        self.params['Title'] = 'Primary flux of ' + tmpparams['particle'] +\
            ' [' + re.sub('#', 'nb particles', tmpparams['flux_unit']) + ']'
        self.params['Xaxis'] = 'Energy[' + tmpparams['energy_unit'] + ']'
        self.params['filename'] = filename
        self.params['interpolation'] = tmpparams['interpolation']
        #set data
        zerocol = zeros(len(tmpdata))
        self.data = column_stack((zerocol, zerocol, tmpdata, zerocol))
        for i in arange(0, len(self.data), 1):
            if i == 0:
                bin_width = self.data[1, 2] - self.data[0, 2]
            else:
                bin_width = self.data[i, 2] - self.data[i - 1, 2]
            self.data[i, 0] = self.data[i, 2] - bin_width / 2
            self.data[i, 1] = self.data[i, 2] + bin_width / 2
        return
        
    def save_as_gps(self, filename, shape = ''):
        """Saves the histogram information in a set of Geant4 GPS gun
        commands that will recreate this spectrum."""
        macrofile = open(filename, 'w')
        macrofile.write('#planetoparse generated GPS source setup for ' +\
            self.particle + '\n')
        macrofile.write('/gps/source/add ' +\
            str(self.params['normalisation_factor']) + '\n')
        if 'primary' in self.particle:
            partname = self.particle.split(' ')[1]
        else:
            partname = self.particle
        macrofile.write('/gps/particle ' + partname + '\n')
        #shape config:
        if not shape == '':
            if shape == 'RAD_above':
                macrofile.write('/gps/pos/type Beam\n')
                macrofile.write('/gps/pos/centre 0. 0. 20. cm\n')
                macrofile.write('/gps/pos/shape Circle\n')
                macrofile.write('/gps/pos/radius 2. cm\n')
                macrofile.write('/gps/direction 0. 0. -1.\n')
            elif shape == 'RAD_below':
                macrofile.write('/gps/pos/type Beam\n')
                macrofile.write('/gps/pos/centre 0. 0. -20. cm\n')
                macrofile.write('/gps/pos/shape Circle\n')
                macrofile.write('/gps/pos/radius 2. cm\n')
                macrofile.write('/gps/direction 0. 0. 1.\n')
            elif shape == 'roveronmars':
                macrofile.write('/gps/pos/type Plane\n')
                macrofile.write('/gps/pos/shape Square\n')
                macrofile.write('/gps/pos/centre 0. 0. -.5 m\n')
                macrofile.write('/gps/pos/halfx .5 m\n')
                macrofile.write('/gps/pos/halfy .5 m\n')
                macrofile.write('/gps/ang/type cos\n')
            else:
                macrofile.write(str(shape) + '\n')
        #histogram points:
        macrofile.write('/gps/ene/type User\n')
        macrofile.write('/gps/hist/type energy\n')
        #first point:
        macrofile.write('/gps/hist/point ' + str(self.data[0, 0]) + '\n')
        for i in arange(0, len(self.data)):
            macrofile.write('/gps/hist/point ' + str(self.data[i, 1]) +\
                ' ' + str(self.data[i, 3]) + '\n')
        macrofile.close()
        return
            
    def save_as_2d_gps(self, filename, shape = '', thetalimits = None,
                       thetaoffset = 0., thetamultiplier = 1.):
        """Saves the histogram information in a set of Geant4 GPS gun
        commands that will recreate this spectrum. Only really works
        with Ekin vs cos theta histograms."""
        macrofile = open(filename, 'w')
        macrofile.write('#planetoparse generated 2D GPS source setup for ' +\
        self.particle + '\n\n')
        #shape config:
        if not shape == '':
            if shape == 'roveronmars':
                shape = ''
                shape += '/gps/pos/type Plane\n'
                shape += '/gps/pos/shape Square\n'
                shape += '/gps/pos/centre 0. 0. -.5 m\n'
                shape += '/gps/pos/halfx .5 m\n'
                shape += '/gps/pos/halfy .5 m\n'
        #get number of bins for energy and cosine theta:
        xedges = []
        for line in self.data:
            if not line[0] in xedges:
                xedges.append(line[0])
            if not line[1] in xedges:
                xedges.append(line[1])
        num_ebins = len(xedges) - 1
        #compute source intensity:
        intensity = self.params['normalisation_factor']
        mask = (self.data[:, 2] < 0.) * (self.data[:, 3] < 0.)        
        total_flux = sum(self.data[:, 4][mask])
        #write source definition
        for i in arange(0, len(self.data), num_ebins):
            data = self.data[i:i + num_ebins]
            mask = (data[:, 2] < 0.) * (data[:, 3] < 0.)
            flux = sum(data[:, 4][mask])
            source_intensity = intensity * flux / total_flux
            if not source_intensity == 0.:
                string = self.__gen_1d_gps(data, shape, source_intensity, 
                                           thetalimits = thetalimits,
                                           thetaoffset = thetaoffset,
                                           thetamultiplier = thetamultiplier)
                macrofile.write(string)
                macrofile.write('\n')
        macrofile.close()
        return
        
    def __gen_1d_gps(self, hist, shape, intensity, thetalimits, thetaoffset,
                     thetamultiplier):
        """Generates a representation of the histogram data as a 1D Geant4
        GPS gun command sequence."""
        #make data copy
        data = empty_like(hist)
        data[:] = hist
        res = ''
        res += '/gps/source/add ' + str(intensity) + '\n'
        res += '/gps/particle ' + self.particle + '\n'
        res += shape
        #histogram points:
        res += '/gps/ene/type User\n'
        res += '/gps/hist/type energy\n'
        elow = data[0, 0]
        ehigh = data[0, 1]
        res += '/gps/hist/point ' + str(elow) + '\n'
        res += '/gps/hist/point ' + str(ehigh) + ' 1.\n'
        res += '/gps/ang/type user\n'
        res += '/gps/hist/type phi\n'
        res += '/gps/hist/point 0.\n'
        res += '/gps/hist/point ' + str(2*pi) + ' 1.\n'
        res += '/gps/hist/type theta\n'
        #prepare the data array:
        mask = ones(len(data[:,0]), dtype = bool)
        if not thetalimits is None:
            thetalimits = cos(array(thetalimits))
            costhetamin = amin(thetalimits)
            costhetamax = amax(thetalimits)
            thetamasklower = (data[:, 2] > costhetamin) * (data[:, 2] < costhetamax)
            thetamaskupper = (data[:, 3] > costhetamin) * (data[:, 3] < costhetamax)
            mask = mask * thetamasklower * thetamaskupper
        data = data[mask]
        data[:, 2] = arccos(data[:, 2]*thetamultiplier + thetaoffset)
        data[:, 3] = arccos(data[:, 3]*thetamultiplier + thetaoffset)
        #If upper angle bin < lower angle bin, swap the two columns:
        if (data[:,3] < data[:,2]).any():
            data = data.T[[0, 1, 3, 2, 4, 5]].T
        data = data[data[:, 3].argsort()]
        if not (data[:, 4] == 0.).all():
            #first point:
            res += '/gps/hist/point ' + str(data[0, 2]) + '\n'
            for i in arange(0, len(data)):
                res += '/gps/hist/point ' + str(data[i, 3]) + ' ' +\
                    str(data[i, 4]) + '\n'
            return res
        else:
            return ''

####################
#planetoparse class definition
####################

class planetoparse:
    """Parses Planetocosmics ASCII output for interactive use. Initialize 
    with filename to parse, see members for parse results. Save and load 
    saves and loads the data to and from a file."""
    
    def __init__(self, filename = None, verbosity = 0, sort_config = None):
        self.primaries = 0
        self.normalisation = ''
        self.params = {}
        self.hists2d = []
        self.cosmonuc = None
        self.hists1d = []
        self.flux_up = {}
        self.flux_down = {}
        self.flux_angular = {}
        self.edep_soil = []
        self.edep_atmo = []
        self.primhists = {}
        self.flux2d_up = {}
        self.flux2d_down = {}
        if sort_config is None:
            try:
                from planetoparse_cfg import sort_config
            except ImportError:
                sort_config = DEFAULT_SORT_CONFIG
        sort_conf_keys = sort_config.keys()
        sort_conf_keys.sort()
        req_keys = DEFAULT_SORT_CONFIG.keys()
        req_keys.sort()
        if not sort_conf_keys == req_keys:
            self.sort_config = DEFAULT_SORT_CONFIG
        else:
            self.sort_config = sort_config
        if not filename is None:
            self.parse_file(filename, verbosity = verbosity)
        return

    def __parse_params(self, line):
        """Parses parameter information. Returns name, value pair if found,
        None, None pair if not found."""
        parsed = re.match(PARAMS_PATTERN, line)
        if parsed:
            name = parsed.group(1)
            try:
                value = float64(parsed.group(2))
            except ValueError:
                value = parsed.group(2)
            return name, value
        else:
            return None, None
        
    def parse_file(self, filename, verbosity = 0):
        """Parse a Planetocosmics ASCII output file."""
        import os
        import sys
        import cStringIO
        import subprocess
        if filename.endswith('.ascii.tar.gz'):    
            p = subprocess.Popen(["tar", "-xOzf", filename], stdout = subprocess.PIPE)
            infile = cStringIO.StringIO(p.communicate()[0])
            assert p.returncode == 0
        elif filename.endswith('.ascii.gz'):          
            p = subprocess.Popen(["zcat", filename], stdout = subprocess.PIPE)
            infile = cStringIO.StringIO(p.communicate()[0])
            assert p.returncode == 0
        elif filename.endswith('.ascii'): 
            infile = open(filename, 'r')
        else:
            print 'ERROR: Unknown file type. Please load .ascii, .tar.gz or .gz files. Continue at own risk'
            infile = open(filename, 'r')
        line = infile.readline()
        while not line == '':
            #parse global information
            while not line == DELIMITER:
                name, value = self.__parse_params(line)
                if not name is None:
                    if name == 'nb_of_primaries':
                        self.primaries = int64(value)
                    elif name == 'normalisation_type':
                        self.normalisation = value
                    else:
                        self.params[name] = value
                line = infile.readline()
            #we have a delimiter, determine if it's a 1d or 2d histogram 
            #and load accordingly:
            line = infile.readline()
            if line.split('\t')[0] == 'Histogram2D':
                infile, line = self.__parse_2d_hist(infile, line)
            else:
                infile, line = self.__parse_1d_hist(infile, line)
        infile.close()
        #print additional output, depending on verbosity
        if verbosity > 0:
            print 'Successfully parsed file' + filename
        if verbosity > 1:
            self.print_stats()
        if verbosity > 2:
            self.print_empty()
        if verbosity > 1:
            print
        return
            
    def print_stats(self):
        """Print information on the number of histograms available."""
        print 'Number of primaries:', self.primaries
        print 'Normalisation:', self.normalisation
        count = 0
        if not self.cosmonuc is None:
            print 'Cosmogenic nuclide histogram present'
            count += 1    
        tmpcount = 0
        for particle in self.flux2d_up:
            for detector in self.flux2d_up[particle]:
                tmpcount += len(self.flux2d_up[particle][detector])
        print 'Upwards 2D flux histograms:', tmpcount
        count += tmpcount
        tmpcount = 0
        for particle in self.flux2d_down:
            for detector in self.flux2d_down[particle]:
                tmpcount += len(self.flux2d_down[particle][detector])
        print 'Downwards 2D flux histograms:', tmpcount
        count += tmpcount
        print 'Other 2D histograms:', len(self.hists2d)
        count += len(self.hists2d)
        tmpcount = 0
        for particle in self.primhists:
            tmpcount += len(self.primhists[particle])
        print 'Primaries histograms:', tmpcount
        count += tmpcount
        print 'Atmosphere energy deposition histograms:', len(self.edep_atmo)
        count += len(self.edep_atmo)
        print 'Soil energy deposition histograms:', len(self.edep_soil)
        count += len(self.edep_soil)
        if not len(self.flux_up) == 0:
            print 'Upward flux histograms:', len(self.flux_up)*\
                len(self.flux_up[self.flux_up.keys()[0]])
            count += len(self.flux_up)*\
                len(self.flux_up[self.flux_up.keys()[0]])
        else:
            print 'No upward flux histograms'
        if not len(self.flux_down) == 0:
            print 'Downward flux histograms:', len(self.flux_down)*\
                len(self.flux_down[self.flux_down.keys()[0]])
            count += len(self.flux_down)*\
                len(self.flux_down[self.flux_down.keys()[0]])
        else:
            print 'No downward flux histograms'
        if not len(self.flux_angular) == 0:
            print 'Angular flux histograms:', len(self.flux_angular)*\
                len(self.flux_angular[self.flux_angular.keys()[0]])
            count += len(self.flux_angular)*\
                len(self.flux_angular[self.flux_angular.keys()[0]])
        else:
            print 'No angular flux histograms'
        print 'Other 1D histograms:', len(self.hists1d)
        count += len(self.hists1d)
        print
        print 'Total:', count, 'histograms'
        return

    def print_empty(self):
        """Print information on empty histograms, if any."""
        def parse_title(hist):
            """Parses and returns title of histogram."""
            if 'Title' in hist.params:
                titleparse = re.match('(.*)\s*\[(.*)\]',
                                      hist.params['Title'])
                res = titleparse.group(1)
            else:
                res = 'Unknown histogram title'
            return res
        message = ''
        count = 0
        #cosmonuc:
        if self.cosmonuc.isempty():
            message += '\tCosmogenic nuclides histogram\n'
            count += 1
        #flux2d_down:
        for particle in self.flux2d_down:
            for detector in self.flux2d_down[particle]:
                for hist in self.flux2d_down[particle][detector]:
                    if hist.isempty():
                        message += '\tflux2d_down, particle ' + particle +\
                            ', detector ' + str(detector) + ':' +\
                            parse_title(hist) + '\n'
                        count += 1
        #flux2d_up:
        for particle in self.flux2d_up:
            for detector in self.flux2d_up[particle]:
                for hist in self.flux2d_up[particle][detector]:
                    if hist.isempty():
                        message += '\tflux2d_up, particle ' + particle +\
                            ', detector ' + str(detector) + ':' +\
                            parse_title(hist) + '\n'
                        count += 1
        #hists2d:
        for hist in self.hists2d:
            if hist.isempty():
                message += '\thists2d: ' + parse_title(hist) + '\n'
                count += 1
        #edep_atmo:
        for hist in self.edep_atmo:
            if hist.isempty():
                message += '\tedep_atmo: ' + parse_title(hist) + '\n'
                count += 1
        #edep_soil:
        for hist in self.edep_soil:
            if hist.isempty():
                message += '\tedep_soil: ' + parse_title(hist) + '\n'
                count += 1
        #primaries:
        for particle in self.primhists:
            for hist in self.primhists[particle]:
                if hist.isempty():
                    message += '\tprimhists, particle ' + particle + '\n'
                    count += 1
        #flux_down:
        for particle in self.flux_down:
            for detector in self.flux_down[particle]:
                if self.flux_down[particle][detector].isempty():
                    message += '\tflux_down, particle ' + particle +\
                        ', detector ' + str(detector) + '\n'
                    count += 1
        #flux_up:
        for particle in self.flux_up:
            for detector in self.flux_up[particle]:
                if self.flux_down[particle][detector].isempty():
                    message += '\tflux_up, particle ' + particle +\
                        ', detector ' + str(detector) + '\n'
                    count += 1
        #flux_angular:
        for particle in self.flux_angular:
            for detector in self.flux_angular[particle]:
                if self.flux_angular[particle][detector].isempty():
                    message += '\tflux_angular, particle ' + particle +\
                        ', detector ' + str(detector) + '\n'
                    count += 1
        #hists1d:
        for hist in self.hists2d:
            if hist.isempty():
                message += '\thists1d: ' + parse_title(hist) + '\n'
                count += 1
        #finalize and print message:
        if message == '':
            message = 'No all-zero histograms detected.'
        else:
            message = 'The following all-zero histograms \
                have been detected:\n' + message
            message += '\nTotal count: ' + str(count)
        print message
        return

    def __parse_hist(self, infile, line):
        """Parses histogram information into a histdata instance."""
        res = histdata()
        line = line.split('\t')
        res.type = line[0]
        res.title = line[1][:-1]
        infile.readline()
        line = infile.readline()
        #parse histogram information
        while not line == DELIMITER:
            name, value = self.__parse_params(line)
            if not name is None:
                res.params[name] = value
            line = infile.readline()
            #check if we have a data line, if so, break the loop
            try:
                dobreak = True
                array(line.split(), dtype = float64)
            except ValueError:
                dobreak = False
            if dobreak:
                break
        tmpdat = []
        while not (line == DELIMITER or line == ''):
            tmpdat.append(array(line.split(), dtype = float64))
            line = infile.readline()
        res.data = array(tmpdat)
        return res, infile, line        

    def __parse_2d_hist(self, infile, line):
        """Parses 2D histogram information into a histdata instance."""
        res, infile, line = self.__parse_hist(infile, line)
        title = res.title.split('/')
        if title[1] == 'COSMONUC':
            self.cosmonuc = res
        elif title[1] == 'FLUX':
            res.particle = title[3]
            res.detector = int(title[2][3:])
            if title[4] in self.sort_config['flux2d_down']:
                if not res.particle in self.flux2d_down:
                    self.flux2d_down[res.particle] = {}
                if not res.detector in self.flux2d_down[res.particle]:
                    self.flux2d_down[res.particle][res.detector] = []
                self.flux2d_down[res.particle][res.detector].append(res)
            elif title[4] in self.sort_config['flux2d_up']:
                if not res.particle in self.flux2d_up:
                    self.flux2d_up[res.particle] = {}
                if not res.detector in self.flux2d_up[res.particle]:
                    self.flux2d_up[res.particle][res.detector] = []
                self.flux2d_up[res.particle][res.detector].append(res)
            else:
                self.hists2d.append(res)
        else:
            self.hists2d.append(res)
        return infile, line
        
    def __parse_1d_hist(self, infile, line):
        """Parses 1D histogram information into a histdata instance."""
        res, infile, line = self.__parse_hist(infile, line)
        title = res.title.split('/')
        if title[1] == 'FLUX':
            res.particle = title[3]
            res.detector = int(title[2][3:])
            if title[4] in self.sort_config['flux_down']:
                if not res.particle in self.flux_down:
                    self.flux_down[res.particle] = {}
                self.flux_down[res.particle][res.detector] = res
            elif title[4] in self.sort_config['flux_up']:
                if not res.particle in self.flux_up:
                    self.flux_up[res.particle] = {}
                self.flux_up[res.particle][res.detector] = res
            else:
                self.hists1d.append(res)
        elif title[1] == 'EDEP':
            self.edep_atmo.append(res)
        elif title[1] == 'SOIL_EDEP':
            self.edep_soil.append(res)
        elif title[1] == 'PRIMARY':
            particle = res.title.split('/')[2]
            res.particle = 'primary ' + particle
            res.detector = 0
            #fix primary particle flux being too low:
            res.data[:, 3:] *= 2.
            if not particle in self.primhists:
                self.primhists[particle] = []
            self.primhists[particle].append(res)
        else:
            self.hists1d.append(res)
        return infile, line
        
    def save(self, filename, gzip_flag = True):
        import gzip
        """Save the contained Planetocosmics result information into a 
        binary gziped file for later use."""
        if gzip_flag: 
            outfile = gzip.open(filename + '.gz', 'wb')
        else:
            outfile = open(filename, 'wb')
        cPickle.dump(self.primaries, outfile)
        cPickle.dump(self.normalisation, outfile)
        cPickle.dump(self.params, outfile)
        cPickle.dump(self.hists2d, outfile)
        cPickle.dump(self.flux2d_up, outfile)
        cPickle.dump(self.flux2d_down, outfile)
        cPickle.dump(self.cosmonuc, outfile)
        cPickle.dump(self.hists1d, outfile)
        cPickle.dump(self.flux_up, outfile)
        cPickle.dump(self.flux_down, outfile)
        cPickle.dump(self.edep_soil, outfile)
        cPickle.dump(self.edep_atmo, outfile)
        cPickle.dump(self.primhists, outfile)
        cPickle.dump(self.flux_angular, outfile)
        outfile.close()
        return
        
    def load(self, filename, print_stats = False):
        """Load the (gziped) Planetocosmics result information from a binary file."""
        import gzip
        if filename.endswith('.gz') and not filename.endswith('.tar.gz'):
            infile = gzip.open(filename, 'rb')
        else:
            infile = open(filename, 'rb')
        self.primaries = cPickle.load(infile)
        self.normalisation = cPickle.load(infile)
        self.params = cPickle.load(infile)
        self.hists2d = cPickle.load(infile)
        self.flux2d_up = cPickle.load(infile)
        self.flux2d_down = cPickle.load(infile)
        self.cosmonuc = cPickle.load(infile)
        self.hists1d = cPickle.load(infile)
        self.flux_up = cPickle.load(infile)
        self.flux_down = cPickle.load(infile)
        self.edep_soil = cPickle.load(infile)
        self.edep_atmo = cPickle.load(infile)
        self.primhists = cPickle.load(infile)
        try:
            self.flux_angular = cPickle.load(infile)
        except EOFError:
            self.flux_angular = {}
        infile.close()
        if print_stats:
            self.print_stats()
        return
        
    def __save_hist_to_ascii(self, hist, outfile):
        """Saves histogram information to an ASCII file"""
        outfile.write(DELIMITER)
        outfile.write(hist.type + '\t' + hist.title + '\n')
        outfile.write(DELIMITER)
        for param in hist.params:
            outfile.write(param + ' : ' + str(hist.params[param]) + '\n')
        savetxt(outfile, hist.data)
        return        
        
    def save_ascii(self, filename, gzip_flag = False):
        """Save the contained Planetocosmics result information into an 
        (gziped) ASCII file for later use."""
        import gzip
        if gzip_flag: 
            outfile = gzip.open(filename + '.ascii.gz', 'wb')
        else:
            outfile = open(filename + '.ascii', 'w')
        if len(self.flux_angular) > 0:
            print 'WARNING: flux_angular will not be saved.'
        outfile.write('nb_of_primaries : ' + str(self.primaries) + '\n')
        outfile.write('normalisation_type : ' + str(self.normalisation) +\
            '\n')
        for param in self.params:
            outfile.write(param + ' : ' + str(self.params[param]) + '\n')
        self.__save_hist_to_ascii(self.cosmonuc, outfile)
        for particle in self.flux2d_up:
            for detector in self.flux2d_up[particle]:
                for hist in self.flux2d_up[particle][detector]:
                    self.__save_hist_to_ascii(hist, outfile)
        for particle in self.flux2d_down:
            for detector in self.flux2d_down[particle]:
                for hist in self.flux2d_down[particle][detector]:
                    self.__save_hist_to_ascii(hist, outfile)
        for hist in self.hists2d:
            self.__save_hist_to_ascii(hist, outfile)
        for particle in self.flux_up:
            for detector in self.flux_up[particle]:
                self.__save_hist_to_ascii(self.flux_up[particle][detector],
                                          outfile)
        for particle in self.flux_down:
            for detector in self.flux_down[particle]:
                self.__save_hist_to_ascii(self.flux_down[particle][detector],
                                          outfile)
        for hist in self.edep_soil:
            self.__save_hist_to_ascii(hist, outfile)
        for hist in self.edep_atmo:
            self.__save_hist_to_ascii(hist, outfile)
        for particle in self.primhists:
            for hist in self.primhists[particle]:
                self.__save_hist_to_ascii(hist, outfile)
        for hist in self.hists1d:
            self.__save_hist_to_ascii(hist, outfile)
        outfile.close()
        return

    def set_scale_per_nuc(self, scale, particle, weight = None):
        """Sets (scale = True) or unsets (scale = False) per nucleus scaling
        for all flux histograms of a given particle."""
        count = 0
        if scale and weight is None:
            print 'ERROR: Need particle weight for scaling.'
            return
        if particle in self.primhists:
            for hist in self.primhists[particle]:
                if scale:
                    hist.scale_per_nuc(weight)
                else:
                    hist.unscale_per_nuc()
                count += 1
        if particle in self.flux_up:
            for detector in self.flux_up[particle]:
                if scale:
                    self.flux_up[particle][detector].scale_per_nuc(weight)
                else:
                    self.flux_up[particle][detector].unscale_per_nuc()
                count += 1
        if particle in self.flux_down:
            for detector in self.flux_down[particle]:
                if scale:
                    self.flux_down[particle][detector].scale_per_nuc(weight)
                else:
                    self.flux_down[particle][detector].unscale_per_nuc()
                count += 1
        if scale:
            print 'Scaled', count, particle, 'flux histograms with weight', \
                weight
        else:
            print 'Unscaled', count, particle, 'flux histograms'
        return
            
        
    def set_scale_per_sterad(self, scale):
        """Sets (scale = True) or unsets (scale = False) per steradian 
        scaling for all flux histograms."""
        count = 0
        for particle in self.primhists:
            for hist in self.primhists[particle]:
                if scale:
                    hist.scale_per_sterad()
                else:
                    hist.unscale_per_sterad()
                count += 1
        for particle in self.flux_up:
            for detector in self.flux_up[particle]:
                if scale:
                    self.flux_up[particle][detector].scale_per_sterad()
                else:
                    self.flux_up[particle][detector].unscale_per_sterad()
                count += 1
        for particle in self.flux_down:
            for detector in self.flux_down[particle]:
                if scale:
                    self.flux_down[particle][detector].scale_per_sterad()
                else:
                    self.flux_down[particle][detector].unscale_per_sterad()
                count += 1
        if scale:
            print 'Scaled', count, 'flux histograms'
        else:
            print 'Unscaled', count, 'flux histograms'
        return

    def set_angular_flux(self, limits=(-1., 1.), detector_levels=[]):
        """
        Sets angular limited flux histos to self.flux_angular.
        pass angular limits tuple with limits=(-1, 1).
        If detector_levels list is empty every detector level will be used (default) or
        if detector_levels=[1,2,3] detector levels 1, 2 and 3 are used for angular conversion.
        """
        from planetotools import project_data
        for hist2d in self.hists2d:
            title = hist2d.title
            if '/DET' in title and '/CosZenVsEkin' in title:
                detector_level = int(re.match('([a-zA-Z]*)([0-9]*)', title.split('/')[2]).group(2))
                if detector_level in detector_levels or len(detector_levels) == 0:
                    print 'Processing histogram:', title
                    element = title.split('/')[3]
                    hist = project_data(hist2d, axis='x', xlimits=None, ylimits=limits)
                    if not element in self.flux_angular:
                        self.flux_angular[element] = {}
                    if not detector_level in self.flux_angular[element]:
                        self.flux_angular[element][detector_level] = histdata(copyhist=hist)
        
    def combine_highz_isotopes(self, verbosity = 0):
        """Combines available high-Z histograms in flux_up and flux_down into
        one histogram per element. High-Z histograms are detected by the 
        string '[0.0]' being present in the particle name. Fluxes are scaled
        to energy/nuc prior to combining, the resulting histogram will have 
        the binning and nuclear weight of the middle isotope."""
        #do this in flux_up, flux_down, flux_angular
        if verbosity > 0:
            print 'Combining high-Z histograms in downward flux...'
        self.__combine_highz_flux(self.flux_down, verbosity = verbosity)
        if verbosity > 0:
            print 'Combining high-Z histograms in upward flux...'
        self.__combine_highz_flux(self.flux_up, verbosity=verbosity)
        if len(self.flux_angular) > 0:
            if verbosity > 0:
                print 'Combining high-Z histograms in angular flux...'
            self.__combine_highz_flux(self.flux_angular, verbosity=verbosity)
        return
        
    def __get_highz_element_list(self, flux_list):
        """Returns a dict of elements with their isotopes as entries."""
        #get list of isotopes
        isotopes = []
        for particle in flux_list.keys():
            if '[0.0]' in particle:
                isotopes.append(re.sub('\[0.0\]', '', particle))
        #get list of elements
        elements = {}
        for isotope in isotopes:
            parse_isotope = re.match('([a-zA-Z]*)([0-9]*)', isotope)
            if not parse_isotope is None:
                if not parse_isotope.group(1) in elements:
                    elements[parse_isotope.group(1)] = []
                if not parse_isotope.group(2) in \
                    elements[parse_isotope.group(1)]:
                    elements[parse_isotope.group(1)].append(
                        parse_isotope.group(2))
        return elements

    def __combine_highz_flux(self, flux_list, verbosity = 0):
        """Combines high-Z flux histograms into one for each element."""
        from . import planetotools as pt
        elements = self.__get_highz_element_list(flux_list)
        if verbosity > 0:
            if len(elements) > 0:
                string = 'Found the following high-Z elements:'
                for element in elements:
                    string += '\n\t' + element + ', isotopes: '
                    for isotope in elements[element]:
                        string += isotope + ', '
                    string = string[:-2]
                print string
            else:
                #there are no elements, might as well just return
                print 'No high-Z elements found.'
                return
        #combine histograms for each element
        for element in elements:
            detectors = self.__get_common_detectors(flux_list, element,
                                                    elements[element])
            if verbosity > 1:
                print 'Combining ' + element + ' histograms in ' +\
                    str(len(detectors)) + ' detector levels...'
            for detector in detectors:
                #find the middle isotope, we want to use that binning:
                middle_index = int64(floor(len(elements[element])/2))
                middle_isotope = elements[element][middle_index]
                #get the list of the remaining isotopes:
                isotopes = elements[element][:middle_index] +\
                    elements[element][middle_index + 1:]
                #copy the histogram of the middle isotope,
                #we will use that for the basis
                middle_hist = flux_list[self.__get_particle_name(
                    element, middle_isotope)][detector]
                middle_hist.scale_per_nuc(float64(middle_isotope))
                res = histdata(copyhist = middle_hist)
                res.params['Title'] = re.sub(self.__get_particle_name(
                    element, middle_isotope, regex = True), element, 
                    res.params['Title'])
                res.particle = element
                #move through remaining isotopes and add up the histograms:
                for isotope in isotopes:
                    hist = flux_list[self.__get_particle_name(element, 
                                                              isotope)]\
                                                                  [detector]
                    hist.scale_per_nuc(float64(isotope))
                    res = pt.__combine_single_hists(res, hist)
                #add the result into the flux list:
                if not element in flux_list:
                    flux_list[element] = {}
                flux_list[element][detector] = res
        return
        
    def __get_particle_name(self, element, isotope, regex = False):
        """Returns the Geant4 particle name for a given element and isotope."""
        if not regex:
            return element + isotope + '[0.0]'
        else:
            return element + isotope

    def __get_common_detectors(self, flux_list, element, isotope_list):
        """Returns a list of detectors that are common for a set of
        histograms."""
        res = flux_list[self.__get_particle_name(element, 
                                                 isotope_list[0])].keys()
        for isotope in isotope_list[1:]:
            tmp = flux_list[self.__get_particle_name(element, 
                                                     isotope)].keys()
            for detector in res:
                if not detector in tmp:
                    res.remove(detector)
        return res
                
    def print_nonempty_highz(self):
        """Print a list of high-Z histograms that are not empty."""
        downlist, downcount = self.__get_nonempty_highz_string(
            self.flux_up, '\tUpward flux histograms')
        uplist, upcount = self.__get_nonempty_highz_string(
            self.flux_up, '\tUpward flux histograms')
        if not (downcount == 0 and upcount == 0):
            output = 'The following non-zero high-Z histograms have been \
                detected:\n' + downlist + uplist
            output += '\n\nTotal count: ' + str(downcount + upcount)
        else:
            output = 'No non-zero high-Z histograms have been detected.'
        print output
        
    def __get_nonempty_highz_string(self, flux_list, prefix):
        """Returns list and count of non-empty high-Z flux histograms."""
        res = ""
        count = 0
        elements = self.__get_highz_element_list(self.flux_up)
        for element in elements:
            for detector in flux_list[element]:
                if not flux_list[element][detector].isempty():
                    res += prefix + ', element ' + element + ', detector ' +\
                        str(detector) + '\n'
                    count += 1
        return res[:-1], count

