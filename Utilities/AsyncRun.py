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

 Title: AsyncRun.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2009-05-01
 Description: A simple wrapper class to initiate threaded processes
 SeeAlso: threading
 Constraints: Untested on single-core processors.

 Version: $Rev: 642 $
 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: AsyncRun.py 642 2012-02-21 07:54:04Z nsummons $
"""

import threading
__version__ = '$Id: AsyncRun.py 642 2012-02-21 07:54:04Z nsummons $'
class AsyncRun(threading.Thread):
    """
    A wrapper around a function to set the function running in a
    separate thread.
    Input: function - name of the function to run
           args - dictionary of arguments to pass to the function
    Output: A Thread object that can be initiated using the start()
            method
    Example:  thread = AsyncRun(function, args)
              thread.start()
              # Main program thread continues here
              # You can also wait for the thread to complete before
              # continuing by using the join() method
              thread.join()
    """
    def __init__(self, function, args):
        threading.Thread.__init__(self)
        self.function = function
        self.args = args

    def run(self):
        self.function(**self.args)

