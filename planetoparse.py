#!/usr/bin/env python

from numpy import *
import re
import cPickle

DELIMITER = '//////////////////////////////////////////////////\n'
PARAMS_PATTERN = '(\w*)\s*:\s*(.*)\n'

####################
#histdata class definition
####################

class histdata:
	"""Provides histogram data to plotting routines. Allows scaling and unscaling per nucleus and per sterad, as well as saving and loading to and from Planetocosmics flux data definition files."""
	
	def __init__(self, copyhist = None):
		self.type = None
		self.title = None
		self.params = {}
		self.data = None
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
			self.data[:,:3] /= self.nuc_weight
			self.data[:,3:] *= self.nuc_weight
			titleparse = re.match('(.*)\s*\[(.*)\]', self.params['Xaxis'])
			self.params['Xaxis'] = titleparse.group(1) + '[' +  titleparse.group(2) + '/nuc]'
			self.scaled_per_nuc = True
			return True
		else:
			return False
			
	def unscale_per_nuc(self):
		"""Remove energy scaling per nucleon from histogram."""
		if self.scaled_per_nuc:
			self.data[:,:3] *= self.nuc_weight
			self.data[:,3:] /= self.nuc_weight
			self.params['Xaxis'] = re.sub('/nuc', '', self.params['Xaxis'])
			self.scaled_per_nuc = False
			return True
		else:
			return False
						
	def scale_per_sterad(self):
		"""Scale histogram to Energy/sr."""
		if not self.scaled_per_sterad:
			self.data[:,3] /= 2*pi
			self.data[:,4] /= 2*pi
			titleparse = re.match('(.*)\s*\[(.*)\]', self.params['Xaxis'])
			self.params['Xaxis'] = titleparse.group(1) + '[' +  titleparse.group(2) + '/sr]'
			self.scaled_per_sterad = True
			return True
		else:
			return False
			
	def unscale_per_sterad(self):
		"""Remove energy scaling per sr from histogram."""
		if self.scaled_per_sterad:
			self.data[:,3] *= 2*pi
			self.data[:,4] *= 2*pi
			self.params['Xaxis'] = re.sub('/sr', '', self.params['Xaxis'])
			self.scaled_per_sterad = False
			return True
		else:
			return False
			
	def isempty(self):
		"""Returns true if the histogram is all-zero, false if it is not. None is returned when no valid histogram is loaded."""
		if self.type == 'Histogram1D':
			return self.__histogram_1d_empty()
		elif self.type == 'Histogram2D':
			return self.__histogram_2d_empty()
		else:
			return None
			
	def __histogram_2d_empty(self):
		return (self.data[:,4] == zeros(len(self.data[:,4]))).all()

	def __histogram_1d_empty(self):
		return (self.data[:,3] == zeros(len(self.data[:,3]))).all()

	def save_as_flux(self, filename):
		"""Save the histogram in a format that can be read in Planetocosmics for primary flux definition. Pass the output filename."""
		if self.type == 'Histogram1D':
			def get_unit(title):
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
				outfile.write(str(self.data[i, 2]) + '\t' + str(self.data[i, 3]) + '\n')
			outfile.close()
			return
		else:
			print 'ERROR: Can only save 1D histogram data as flux definition.'
			return
			
	def save_data(self, filename):
		"""Save only the data array to disk."""
		savetxt(filename, self.data)
		return
			
	def __parse_params(self, line):
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
		"""Load the histogram data from a file containing Planetocosmics primary flux definitions."""
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
		self.params['Title'] = 'Primary flux of ' + tmpparams['particle'] + ' [' + re.sub('#', 'nb particles', tmpparams['flux_unit']) + ']'
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
			
####################
#planetoparse class definition
####################

class planetoparse:
	"""Parses Planetocosmics ASCII output for interactive use. Initialize with filename to parse, see members for parse results. Save and load saves and loads the data to and from a file."""
	
	def __init__(self, filename = None, verbosity = 0):
		self.primaries = 0
		self.normalisation = ''
		self.params = {}
		self.hists2d = []
		self.cosmonuc = None
		self.hists1d = []
		self.flux_up = {}
		self.flux_down = {}
		self.edep_soil = []
		self.edep_atmo = []
		self.primhists = {}
		if not filename is None:
			self.parse_file(filename, verbosity = verbosity)
		return
		
	def __parse_params(self, line):
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
		"""Parse a Planetocosmics ASCII output file. Set print_stats to True to get information on the parsing results."""
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
			#we have a delimiter, determine if it's a 1d or 2d histogram and load accordingly:
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
		print 'Other 2D histograms:', len(self.hists2d)
		count += len(self.hists2d)
		print 'Primaries histograms:', len(self.primhists)
		count += len(self.primhists)
		print 'Atmosphere energy deposition histograms:', len(self.edep_atmo)
		count += len(self.edep_atmo)
		print 'Soil energy deposition histograms:', len(self.edep_soil)
		count += len(self.edep_soil)
		if not len(self.flux_up) == 0:
			print 'Upward flux histograms:', len(self.flux_up)*len(self.flux_up[self.flux_up.keys()[0]])
			count += len(self.flux_up)*len(self.flux_up[self.flux_up.keys()[0]])
		else:
			print 'No upward flux histograms'
		if not len(self.flux_down) == 0:
			print 'Downward flux histograms:', len(self.flux_down)*len(self.flux_down[self.flux_down.keys()[0]])
			count += len(self.flux_down)*len(self.flux_down[self.flux_down.keys()[0]])
		else:
			print 'No downward flux histograms'
		print 'Other 1D histograms:', len(self.hists1d)
		count += len(self.hists1d)
		print
		print 'Total:', count, 'histograms'
		return

	def print_empty(self):
		"""Print information on empty histograms, if any."""
		def parse_title(hist):
			if 'Title' in hist.params:
				titleparse = re.match('(.*)\s*\[(.*)\]', hist.params['Title'])
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
			if self.primhists[particle].isempty():
				message += '\tPrimary particle histograms, particle ' + particle + '\n'
				count += 1
		#flux_down:
		for particle in self.flux_down:
			for detector in self.flux_down[particle]:
				if self.flux_down[particle][detector].isempty():
					message += '\tDownward flux histograms, particle ' + particle + ', detector ' + str(detector) + '\n'
					count += 1
		#flux_up:
		for particle in self.flux_up:
			for detector in self.flux_up[particle]:
				if self.flux_down[particle][detector].isempty():
					message += '\tUpward flux histograms, particle ' + particle + ', detector ' + str(detector) + '\n'
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
			message = 'The following all-zero histograms have been detected:\n' + message
			message += '\nTotal count: ' + str(count)
		print message
		return

	def __parse_hist(self, infile, line):
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
		res, infile, line = self.__parse_hist(infile, line)
		if res.title.split('/')[1] == 'COSMONUC':
			self.cosmonuc = res
		else:
			self.hist2d.append(res)
		return infile, line
		
	def __parse_1d_hist(self, infile, line):
		res, infile, line = self.__parse_hist(infile, line)
		title = res.title.split('/')
		if title[1] == 'FLUX':
			res.particle = title[3]
			res.detector = int(title[2][3:])
			if title[4] == '1':
				if not res.particle in self.flux_down:
					self.flux_down[res.particle] = {}
				self.flux_down[res.particle][res.detector] = res
			else:
				if not res.particle in self.flux_up:
					self.flux_up[res.particle] = {}
				self.flux_up[res.particle][res.detector] = res
		elif title[1] == 'EDEP':
			self.edep_atmo.append(res)
		elif title[1] == 'SOIL_EDEP':
			self.edep_soil.append(res)
		elif title[1] == 'PRIMARY':
			particle = res.title.split('/')[2]
			res.particle = 'primary ' + particle
			res.detector = 0
			self.primhists[particle] = res
		else:
			self.hists1d.append(res)
		return infile, line
		
	def save(self, filename):
		"""Save the contained Planetocosmics result information into a binary file for later use."""
		outfile = open(filename, 'wb')
		cPickle.dump(self.primaries, outfile)
		cPickle.dump(self.normalisation, outfile)
		cPickle.dump(self.params, outfile)
		cPickle.dump(self.hists2d, outfile)
		cPickle.dump(self.cosmonuc, outfile)
		cPickle.dump(self.hists1d, outfile)
		cPickle.dump(self.flux_up, outfile)
		cPickle.dump(self.flux_down, outfile)
		cPickle.dump(self.edep_soil, outfile)
		cPickle.dump(self.edep_atmo, outfile)
		cPickle.dump(self.primhists, outfile)
		outfile.close()
		return
		
	def load(self, filename, print_stats = False):
		"""Load the Planetocosmics result information from a binary file."""
		infile = open(filename, 'rb')
		self.primaries = cPickle.load(infile)
		self.normalisation = cPickle.load(infile)
		self.params = cPickle.load(infile)
		self.hists2d = cPickle.load(infile)
		self.cosmonuc = cPickle.load(infile)
		self.hists1d = cPickle.load(infile)
		self.flux_up = cPickle.load(infile)
		self.flux_down = cPickle.load(infile)
		self.edep_soil = cPickle.load(infile)
		self.edep_atmo = cPickle.load(infile)
		self.primhists = cPickle.load(infile)
		infile.close()
		if print_stats:
			self.print_stats()
		return
		
	def __save_hist_to_ascii(self, hist, outfile):
		outfile.write(DELIMITER)
		outfile.write(hist.type + '\t' + hist.title + '\n')
		outfile.write(DELIMITER)
		for param in hist.params:
			outfile.write(param + ' : ' + str(hist.params[param]) + '\n')
		savetxt(outfile, hist.data)
		return		
		
	def save_ascii(self, filename):
		"""Save the contained Planetocosmics result information into an ASCII file for later use."""
		outfile = open(filename, 'w')
		outfile.write('nb_of_primaries : ' + str(self.primaries) + '\n')
		outfile.write('normalisation_type : ' + str(self.normalisation) + '\n')
		for param in self.params:
			outfile.write(param + ' : ' + str(self.params[param]) + '\n')
		self.__save_hist_to_ascii(self.cosmonuc, outfile)
		for hist in self.hists2d:
			self.__save_hist_to_ascii(hist, outfile)
		for particle in self.flux_up:
			for detector in self.flux_up[particle]:
				self.__save_hist_to_ascii(self.flux_up[particle][detector], outfile)
		for particle in self.flux_down:
			for detector in self.flux_down[particle]:
				self.__save_hist_to_ascii(self.flux_down[particle][detector], outfile)
		for hist in self.edep_soil:
			self.__save_hist_to_ascii(hist, outfile)
		for hist in self.edep_atmo:
			self.__save_hist_to_ascii(hist, outfile)
		for particle in self.primhists:
			self.__save_hist_to_ascii(self.primhists[particle], outfile)
		for hist in self.hists1d:
			self.__save_hist_to_ascii(hist, outfile)
		outfile.close()
		return

	def set_scale_per_nuc(self, scale, particle, weight = None):
		"""Sets (scale = True) or unsets (scale = False) per nucleus scaling for all flux histograms of a given particle."""
		count = 0
		if scale and weight is None:
			print 'ERROR: Need particle weight for scaling.'
			return
		if particle in self.primhists:
			if scale:
				self.primhists[particle].scale_per_nuc(weight)
			else:
				self.primhists[particle].unscale_per_nuc()
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
			print 'Scaled', count, particle, 'flux histograms with weight', weight
		else:
			print 'Unscaled', count, particle, 'flux histograms'
		return
			
		
	def set_scale_per_sterad(self, scale):
		"""Sets (scale = True) or unsets (scale = False) per steradian scaling for all flux histograms."""
		count = 0
		for particle in self.primhists:
			if scale:
				self.primhists[particle].scale_per_sterad()
			else:
				self.primhists[particle].unscale_per_sterad()
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

