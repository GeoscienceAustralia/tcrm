#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

 Title: error.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 06/05/08 11:32:AM
 Description:

 Version :$Rev: 512 $

 $Id: error.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, sys, pdb, logging, traceback
filename = os.getenv('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

logger = logging.getLogger()

def errDieWithLog(message=None):
    tb = sys.exc_info()[2]
    stack = []
    while tb:
        stack.append(tb.tb_frame)
        tb = tb.tb_next

    #traceback.print_exc()
    #logger.critical("Locals by frame, innermost last")
    for frame in stack:
        logger.critical("Frame %s in %s at line %s" % (frame.f_code.co_name,
                                                        frame.f_code.co_filename,
                                                        frame.f_lineno))
        for key, value in frame.f_locals.items():
            logger.critical("%s = %s"%(key, repr(value)))

    if message:
        logger.critical(message)
    sys.exit(1)

class errFileOpenError(Exception):
    def __init__(self, fileName):
        self.fileName = fileName
    def __str__(self):
        logger.exception("File open error: cannot open %s"%(repr(self.fileName)))
        return "File open error : cannot open %s"%(repr(self.fileName))
    pass

class errFileCloseError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        logger.exception("File close error: cannot close %s"%(repr(self.fileName)))
        return "File close error: cannot close %s"%(repr(self.fileName))
    pass

class errNetCDFError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        logger.exception("Error in nctools: %s"%repr(self.value))
        return "Error in nctools: %s"%repr(self.value)
    pass



