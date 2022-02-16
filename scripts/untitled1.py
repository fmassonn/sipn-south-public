#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 13:24:06 2022

@author: massonnetf
"""

import numpy as np
import os
from   matplotlib import font_manager



# Load fonts
# Change font globally
# --------------------
#font_dirs  = ["/System/Library/Fonts/", ]
#font_files = font_manager.findSystemFonts(fontpaths = font_dirs)
#font_list  = font_manager.createFontList(font_files)
#font_manager.fontManager.ttflist.extend(font_list)
matplotlib.rcParams['font.family'] = "Arial Narrow"