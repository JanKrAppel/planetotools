#!/usr/bin/env python

from planetotools import *

from os import environ as env
if 'MCD_DATA_DIR' in env:
    from mcdtools import *

