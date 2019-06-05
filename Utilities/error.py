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

 Title: error.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 06/05/08 11:32:AM
 Description:

 Version :$Rev: 642 $

 $Id: error.py 642 2012-02-21 07:54:04Z nsummons $
"""
import sys
import logging

LOG = logging.getLogger()

def errDieWithLog(message=None):
    """
    Capture exception message, log it and die gracefully.

    :param str message: Optional additional error message.

    """

    tb = sys.exc_info()[2]
    stack = []
    while tb:
        stack.append(tb.tb_frame)
        tb = tb.tb_next

    for frame in stack:
        LOG.critical("Frame %s in %s at line %s", frame.f_code.co_name,
                     frame.f_code.co_filename, frame.f_lineno)
        for key, value in list(frame.f_locals.items()):
            LOG.critical("%s = %s", key, repr(value))

    if message:
        LOG.critical(message)
    sys.exit(1)

class ErrFileOpenError(Exception):
    """
    Handle errors when attempting to open files.

    """
    def __init__(self, fileName):
        Exception.__init__()
        self.fileName = fileName
    def __str__(self):
        LOG.exception("File open error: cannot open %s", repr(self.fileName))
        return "File open error : cannot open %s"%(repr(self.fileName))

class ErrFileCloseError(Exception):
    """
    Handle errors when attempting to close files.

    """
    def __init__(self, value):
        Exception.__init__()
        self.value = value
    def __str__(self):
        LOG.exception("File close error: cannot close %s", repr(self.fileName))
        return "File close error: cannot close %s"%(repr(self.fileName))

class ErrNetCDFError(Exception):
    """
    Handle errors when working with netCDF files.

    """
    def __init__(self, value):
        Exception.__init__()
        self.value = value
    def __str__(self):
        LOG.exception("Error in nctools: %s", repr(self.value))
        return "Error in nctools: %s"%repr(self.value)



