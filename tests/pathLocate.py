import os, sys, pdb, logging
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

Title: pathLocator.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au  (Adapted from code snippet found on online forum)
CreationDate: 18-10-2011
Description: A function to determine the directory one level above the Utilities folder.
             Designed to work even if the code has been compiled with py2exe.  This results in the modules
             being built-in to the interpreter, so the path of the executable is returned instead.
"""

def is_frozen():
    # Modules built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def getRootDirectory():
    encoding = sys.getfilesystemencoding()
    if is_frozen():
        return os.path.dirname(str(sys.executable, encoding))
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

def getUnitTestDirectory():
    if is_frozen():
        return os.path.join(getRootDirectory())
    else:
        return os.path.join(getRootDirectory(), 'tests')
