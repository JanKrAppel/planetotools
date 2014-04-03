#!/usr/bin/env python

"""A general-purpose configuration file parser loader function.

24.10.11 - added save_params
19.10.11 - added support for unsetting parameters
13.10.11 - added support for boolean parameters (no '=' in parameter)

Author: Jan Kristoffer Appel
    Email: appel@physik.uni-kiel.de"""

from numpy import *
import re
from os.path import exists

def load_params(file, global_file = None):
    """Parses configuration information from a file and returns a dict containing the entries.
    
    Lines beginning with '#' are interpreted as comments. Each line can only contain
    one parameter.
    Number parameters are returned as numpy.float64 dict entries, every other parameter
    is returned as string dict entry.
    Parameters not followed by a '=' and a value are interpreted as boolean and set to True.
    An entry starting with 'unset' followed by the name of a previously set parameter, e.g.
    'unset steps', will result in that parameter to be removed from the parameter dict. Use
    this to unset parameters set in a global file.  
    If an error occurs, None is returned.
    If global_file is passed, it is parsed before parsing the 'normal' file, thus, 
    any parameters in the 'normal' file can overwrite global ones."""
    params = {}
    if not global_file == None:
        if exists(global_file):
            params = parse_file(global_file)
    if not exists(file):
        print 'ERROR: load_params: File ' + file + ' not found.'
        return None
    else:
        params.update(parse_file(file))
        return params
    
def save_params(params, filename):
    file = open(filename, 'w')
    for key, value in enumerate(params):
        if type(value) == bool:
            file.write(key + '\n')
        else:
            file.write(key + ' = ' + str(value) + '\n')
    file.close()
    
def parse_file(file):
    params = {}
    file = open(file, 'r')
    for line in file:        
        #remove comments:
        line = re.sub('#.*?$', '', line)
        #the [:-1] strips out the last character, which by definition is \n
        if not line[:-1] == '':
            #look for 'normal' parameters:
            matched = False
            values = re.search('([a-zA-Z0-9_-]+)\s*=\s*(.*)', line)
            if not values is None:
                if not values.group(1) == '' and not values.group(2) == '':
                    try:
                        #if it's convertible to a float, add it as a float:
                        params[values.group(1).lower()] = float64(values.group(2))
                        matched = True
                    except ValueError:
                        #if not, just add the string and let the program cope by itself:
                        params[values.group(1).lower()] = values.group(2)
                        matched = True
                else:
                    #if we got an empty string in one of the groups, something went wrong parsing it.
                    print 'ERROR: load_params: Unable to parse entry \'' + line + '\''
            #look for boolean parameters: 
            values = re.search('^([a-zA-Z0-9_-]+)$', line)
            if not values is None:
                params[values.group(1).lower()] = True
                matched = True
            #look for unset commands:
            values = re.search('^unset\s+(.*?)$', line)
            if not values is None:
                if values.group(1).lower() in params:
                    print 'WARNING: Unset parameter ' + values.group(1).lower()
                    params.pop(values.group(1).lower())
                    matched = True
    file.close()
    return params
