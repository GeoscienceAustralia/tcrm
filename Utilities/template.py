"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Title: template.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 03/27/08 8:05:AM
 Description: Read in lines of a text file and replace a given string
 with another value.  e.g. replace all occurences of {KEY} with 0001

 Version :$Rev: 680 $

 $Id: template.py 680 2012-03-19 01:52:37Z carthur $
"""
import os, sys, pdb, logging

__version__ = '$Id: template.py 680 2012-03-19 01:52:37Z carthur $'

import re

def replace(infile, outfile, replacements):
    """replace(infile, outfile, replacements)
    Replace all instances of the keys with values in infile and write to
    outfile.
    Keywords to be replaced in infile should be written as {keyword}
    """
    fi = open(infile, 'r')
    fo = open(outfile, 'w')
    for line in fi:
        newline = line
        for k,v in replacements.items():
            newline = re.sub('{'+k+'}', v, newline)
            
        fo.write(newline)
    fi.close()
    fo.close()
    return True
