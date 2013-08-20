#!/usr/bin/env python

from fmcd import call_mcd, julian
import params_parser
from os import environ
from numpy import *
from MCDConfig import *
from dateutil.parser import parse as parse_date
import re

####################
#MCDToPlanetocosmics .cfg file parser
####################

class mcdparse:
	"""Provides interactive access to MCD .cfg files for MCDToPlanetocosmics."""
	def __init__(self, configfile = None, verbosity = 0):
		self.__params_loaded = False
		self.__data_loaded = False
		self.params = {}
		self.data = {}
		if not configfile is None:
			self.load_params(configfile)
			if self.__params_loaded:
				self.load_data(verbosity = verbosity)
				return
			else:
				print 'ERROR: Unable to read parameters from', configfile
				self.__data_loaded = False
				return
		return
			
	def load_data(self, verbosity = 0):
		"""Builds the composition profile. Parameters need to be loaded first."""
		if not self.__params_loaded:
			print 'ERROR: No parameters loaded, aborting'
			return
		if verbosity > 0:
			self.print_summary()
		#read data		
		self.data = {}
		dat_xz = []
		dat_temp = []
		dat_pres = []
		dat_dens = []
		dat_xz_comp = {}
		for component in self.params['atmo_components'] + self.params['dust_components']:
			dat_xz_comp[component] = []
		for xz in arange(self.params['max_height']*1e3 + self.params['height_steps']*1e3, self.params['min_height']*1e3, -1*self.params['height_steps']*1e3):
			xz, temp, pres, dens, xz_comp = self.__read_xz_composition(xz)
			dat_xz.append(xz)
			dat_temp.append(temp)
			dat_pres.append(pres)
			dat_dens.append(dens)
			for component in dat_xz_comp:
				dat_xz_comp[component].append(xz_comp[component])
		self.data['xz'] = dat_xz
		self.data['temp'] = dat_temp
		self.data['pres'] = dat_pres
		self.data['dens'] = dat_dens
		self.data['comp'] = dat_xz_comp
		self.__build_shield_depth()
		if verbosity > 0:
			print 'Data successfully loaded.'
		self.__data_loaded = True
		return
		
	def __build_shield_depth(self):
		res = []
		for i in arange(0, len(self.data['xz']), 1):
			if i == 0:
				res.append(float64(0.))
			else:
				dh = (self.data['xz'][i - 1] - self.data['xz'][i])*1e5
				res.append((self.data['dens'][i] * dh) + res[-1])
		self.data['shield_depth'] = array(res, dtype = float64)
		return
			
	def print_summary(self):
		"""Print parameters summary message."""
		if self.__params_loaded:
			print 'Parsed file:', self.params['filename']
			for value in dust_scenarios:
				if dust_scenarios[value] == self.params['dust']:
					dustscen = value
			print 'Dust scenario:\t' + dustscen
			if self.params['datekey'] == 0:
				print 'Julian date:\t' + str(self.params['xdate'])
			else:
				print 'Date:\t' + str(self.params['xdate']) + ' Ls, ' + str(self.params['loct']) + ' local time'
			print 'Latitude:\t' + str(self.params['lat'])
			print 'Longitude:\t' + str(self.params['lon'])
			atmo_comp_string = 'Selected atmospheric components:\t'
			for comp in self.params['atmo_components']:
				atmo_comp_string += comp + ' '
			print atmo_comp_string
			if not self.params['dust_components'] == []:
				dust_comp_string = 'Selected dust components:\t'
				for comp in self.params['dust_components']:
					dust_comp_string += comp + ' '
				print dust_comp_string
			print 'Height range from', self.params['min_height'], 'km to', self.params['max_height'], 'km in steps of', self.params['height_steps'], 'km'
			if self.__data_loaded:
				print 'Data successfully loaded.'
		else:
			print 'No parameters loaded.'
		return
		
	def __read_xz_composition(self, xz):
		"""Reads atmospheric composition at given xz."""
		R = 8.3144621
		xz_comp = {}
		(pres, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr) = call_mcd(self.params['zkey'], xz, self.params['lon'], self.params['lat'], self.params['hrkey'], self.params['datekey'], self.params['xdate'], self.params['loct'], self.params['dset'], self.params['dust'], self.params['perturkey'], self.params['seedin'], self.params['gwlength'], ones(100))
		#add normal atmospheric components:
		for comp in self.params['atmo_components']:
			xz_comp[comp] = extvar[component_indices[comp]]
		#add dust components
		Rs = extvar[52]
		air_mass = pres/(Rs*temp)
		dust_mass_ratio = extvar[37]	
		dust_mass = air_mass*dust_mass_ratio
		dust_vol = dust_mass/self.params['dust_density']
		for comp in self.params['dust_components']:
			xz_comp[comp] = dust_vol*self.params['dust_component_ratios'][comp]
		#set surface height and return values
		self.params['surface_height'] = extvar[3]/1e3
		return (extvar[1]/1e3, temp, pres, dens, xz_comp)
		
	def load_params(self, configfile):
		"""Loads parameters from a config file."""
		self.params = params_parser.load_params(configfile, global_file = MCD_GLOBAL_CONFIG)
		#set database location from various sources:
		#1. config file, if not found, then
		#2. environment variable, if not found, then
		#3. fall back on hardcoded location (MCDConfig.py)
		if not 'dset' in self.params:
			if 'MCD_DATA_DIR' in environ:
				self.params['dset'] = environ['MCD_DATA_DIR']
			else:
				self.params['dset'] = MCD_DATA_DIR
		#check if needed parameters are present
		needed_vars = ['zkey', 'lon', 'lat', 'hrkey', 'datekey', 'xdate', 'dust', 'perturkey', 'seedin', 'gwlength', 'min_height', 'max_height', 'height_steps', 'atmo_components']
		for var in needed_vars:
			if not var in self.params:
				print 'ERROR: Missing variable', var, 'in parameter set.'
				self.__params_loaded = False
				return
		#parse dust scenario, if needed
		if type(self.params['dust']) == str:
			if self.params['dust'] in dust_scenarios:
				self.params['dust'] = dust_scenarios[self.params['dust']]
			else:
				print 'ERROR: Invalid dust scenario selected.'
				self.__params_loaded = False
				return
		#check and set date parameters
		if self.params['datekey'] == 0:
			parsed_date = parse_date(self.params['xdate'])
			(ierr, xdate) = julian(parsed_date.month, parsed_date.day, parsed_date.year, parsed_date.hour, parsed_date.minute, parsed_date.second)
			self.params['xdate'] = xdate
			self.params['loct'] = 0
		else:
			if not 'loct' in params:
				print 'ERROR: Argument loct needed when martian time is selected.'
				self.__params_loaded = False
				return
		#set atmospheric composition
		new_atmo_components = []
		for component in re.split('([a-zA-Z0-9]*)', self.params['atmo_components']):
			if not component in ['', ' ']:
				if component in component_indices:
					new_atmo_components.append(component)
		self.params['atmo_components'] = new_atmo_components
		#set dust composition if wanted, else, set to empty
		if 'dust_components' in self.params:
			new_dust_components = []
			new_dust_component_ratios = {}
			for component in re.split('([a-zA-Z0-9]*:[0-9.]*)', self.params['dust_components']):		
				if not component in ['', ' ']:
					comp, ratio = component.split(':')
					new_dust_components.append(comp)
					new_dust_component_ratios[comp] = float64(ratio)
			self.params['dust_components'] = new_dust_components	
			self.params['dust_component_ratios'] = new_dust_component_ratios
			if not 'dust_density' in params:
				print 'WARNING: No dust density specified, defaulting to 1.1 kg/m^3'
				self.params['dust_density'] = 1.1
		else:
			self.params['dust_components'] = []
			self.params['dust_component_ratios'] = {}
			self.params['dust_density'] = 1.
		self.params['filename'] = configfile
		self.__params_loaded = True
		return
		
	def save_params(self, outputfile):
		"""Save parameters into a config file."""
		save_params(self.params, outputfile)
		return
		
	def get_coords_scenario_str(self):
		res = ''
		res += str(self.params['lat']) + ' N, '
		res += str(self.params['lon']) + ' E, '
		for value in dust_scenarios:
			if dust_scenarios[value] == self.params['dust']:
				dustscen = value
		res += dustscen
		return res
		
####################
#Atmotable file parser
####################

COMMENTS_PATTERN = '\s*(\w[a-zA-Z_ ]*\w)\s*:\s*(.*)\n'
DEFINITION_PATTERN = '\s*\\\\(\w*)\{(.*)\}\n'

class atmoparse:
	"""Provides interactive access to Planetocosmics Atmotable files."""
	def __init__(self, atmofile = None, verbosity = 0):
		self.__data_loaded = False
		self.params = {}
		self.data = {}
		if not atmofile is None:
			self.load_data(atmofile = atmofile, verbosity = verbosity)
			return
		return
		
	def __build_shield_depth(self):
		res = []
		for i in arange(0, len(self.data['xz']), 1):
			if i == 0:
				res.append(float64(0.))
			else:
				dh = (self.data['xz'][i - 1] - self.data['xz'][i])*1e5
				res.append((self.data['dens'][i] * dh) + res[-1])
		self.data['shield_depth'] = array(res, dtype = float64)
		return
			
	def load_data(self, atmofile, verbosity = 0):
		"""Builds the composition profile."""
		if self.__data_loaded:
			self.data = {}
			self.params = {}
			self.__data_loaded = False
		self.params['filename'] = atmofile
		self.params['min_height'] = float64(0.)
		dat_xz = []
		dat_temp = []
		dat_pres = []
		dat_dens = []
		dat_xz_comp = {}
		infile = open(atmofile, 'r')
		line = infile.readline()
		####################
		#read comments section
		####################
		while not '\\comments' in line:
			line = infile.readline()
		#read the first comments line
		line = infile.readline()
		while not '\\definition' in line:
			name, value = self.__parse_params(line, COMMENTS_PATTERN)
			if name is None:
				if not 'title' in self.params:
					self.params['title'] = line[1:-1]
			else:
				#translate similar to MCD data
				if name == 'longitude':
					name = 'lon'
				if name == 'latitude':
					name = 'lat'
				if name == 'ground altitude':
					name = 'surface_height'
				self.params[name] = value
			line = infile.readline()
		####################
		#read definition section
		####################
		line = infile.readline()
		while not '\\data' in line:
			name, value = self.__parse_params(line, DEFINITION_PATTERN)
			if not name is None:
				if name == 'ground_altitude':
					name = 'surface_height'
				if name == 'top_altitude':
					name = 'max_height'
				self.params[name] = value
			line = infile.readline()
		####################
		#read data section
		####################
		line = infile.readline()
		#get column mapping
		headers = line.split()
		value_map = {}
		atmo_components = []
		for i in arange(0, len(headers), 1):
			if headers[i] == 'altitude':
				value_map['xz'] = i
			elif headers[i] == 'temperature':
				value_map['temp'] = i
			elif headers[i] == 'pressure':
				value_map['pres'] = i
			elif headers[i] == 'density':
				value_map['dens'] = i
			else:
				atmo_components.append(headers[i])
				value_map[headers[i]] = i
		self.params['atmo_components'] = atmo_components
		for component in atmo_components:
			dat_xz_comp[component] = []
		#read data lines
		line = infile.readline()[:-1]
		while not line == '':
			values = array(line.split('\t')[1:], dtype = float64)
			dat_xz.append(values[value_map['xz']])
			dat_temp.append(values[value_map['temp']])
			dat_pres.append(values[value_map['pres']])
			dat_dens.append(values[value_map['dens']])
			for component in atmo_components:
				dat_xz_comp[component].append(values[value_map[component]])
			line = infile.readline()[:-1]
		#close data file
		infile.close()		
		#set data dict		
		self.data['xz'] = dat_xz
		self.data['temp'] = dat_temp
		self.data['pres'] = dat_pres
		self.data['dens'] = dat_dens
		self.data['comp'] = dat_xz_comp
		self.__build_shield_depth()
		self.__data_loaded = True
		if verbosity > 0:
			self.print_summary()
			print 'Data successfully loaded.'
		return
						
	def __parse_params(self, line, pattern):
		parsed = re.match(pattern, line)
		if parsed:
			name = parsed.group(1)
			try:
				value = float64(parsed.group(2))
			except ValueError:
				value = parsed.group(2)
			return name, value
		else:
			return None, None

	def get_coords_scenario_str(self):
		res = ''
		res += self.params['title'] + ', '
		res += str(self.params['lat']) + ' N, '
		res += str(self.params['lon']) + ' E'
		return res		
		
	def print_summary(self):
		"""Print data summary message."""
		if self.__data_loaded:
			print 'Parsed file:', self.params['filename']
			if 'julian date' in self.params:
				print 'Julian date:\t' + str(self.params['julian date'])
			elif 'xdate' in self.params:
				print 'Date:\t' + str(self.params['xdate']) + ' Ls, ' + str(self.params['loct']) + ' local time'
			elif 'date' in self.params:
				print 'Date:\t' + str(self.params['date'])
			print 'Latitude:\t' + str(self.params['lat'])
			print 'Longitude:\t' + str(self.params['lon'])
			atmo_comp_string = 'Selected atmospheric components:\t'
			for comp in self.params['atmo_components']:
				atmo_comp_string += comp + ' '
			print atmo_comp_string
			print 'Height range from', self.params['min_height'], 'km to', self.params['max_height'], 'km'
		else:
			print 'No data loaded.'
		return
