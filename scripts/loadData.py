#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 13:05:16 2022

@author: massonnetf
"""

from datetime import datetime
from datetime import timedelta


import numpy as np

# Script to load the SIPN south data - forecasts and verification
# The years start in 2015 since this is the first year available
# for observations in the two operational near-real time products
# ---------------------------------------------------------------
exec(open("./namelist_contributors.py").read())

iniDates  = namelist[0]
nIniDates = len(iniDates)



# Load the namelist file
# ----------------------





