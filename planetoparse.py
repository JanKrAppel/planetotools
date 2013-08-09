#!/usr/bin/env python

from numpy import *
import re

DELIMITER = '//////////////////////////////////////////////////\n'
PARAMS_PATTERN = '(\w*)\s*:\s*(.*)\n'

class histdata:
	"""Dummy class to provide histogram encapsulation"""
	def __init__(self, copyhist = None):
		self.title = None
		self.params = {}
		self.data = None
		self.particle = ''
		self.detector = -1
		if not copyhist is None:
			self.title = copyhist.title
			self.params = copyhist.params
			self.data = copyhist.data
			self.particle = copyhist.particle
			self.detector = copyhist.detector
			
class planetoparse:
	"""Parses Planetocosmics ASCII output for interactive use. Initialize with filename to parse, see members for parse results."""
	
	def __init__(self, filename = None, print_stats = False):
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
			self.parse_file(filename, print_stats)
		
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
		
	def parse_file(self, filename, print_stats = False):
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
		if print_stats:
			print 'Successfully parsed', filename
			print 'Parsing summary:'
			print 'Number of primaries:', self.primaries
			print 'Normalisation:', self.normalisation
			if not self.cosmonuc is None:
				print 'Cosmogenic nuclide histogram present'
			print 'Other 2D histograms:', len(self.hists2d)
			print 'Primaries histograms:', len(self.primhists)
			print 'Atmosphere energy deposition histograms:', len(self.edep_atmo)
			print 'Soil energy deposition histograms:', len(self.edep_soil)
			print 'Flux histograms:', len(self.flux_up)*len(self.flux_up[self.flux_up.keys()[0]]), 'up,', len(self.flux_down)*len(self.flux_down[self.flux_up.keys()[0]]), 'down'
			print 'Other 1D histograms:', len(self.hists1d)

	def __parse_hist(self, infile, line):
		res = histdata()
		res.title = line.split('\t')[1][:-1]
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
		while not line == DELIMITER and not line == '':
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
			
if __name__ == '__main__':
	from sys import argv
	results = planetoparse(argv[1], print_stats = True)
