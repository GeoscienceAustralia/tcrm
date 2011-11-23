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

Title: progressbar.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2011-10-14
Description: Draws progress bar on terminal window
""" 
import os, sys, pdb, logging
import time

class ProgressBar():
    def __init__(self, modname):
        self.modname = modname + " "
        self.lastPercentage = None
        self.screenWidth = 79
        self.barWidth = self.screenWidth - len(self.modname) - len(self._getTimeStr()) - 8
        self.start_time = time.time()
        self.update(0.0)

    def update(self, progress, startPos=0, endPos=1):
        prg = progress * (endPos - startPos) + startPos
        if self._percentage(prg) != self.lastPercentage:
            barfill = round(self.barWidth * prg)
            barString = ''.join(['#' for i in range(barfill)] + [' ' for i in range(self.barWidth - barfill)])        
            self.secondsElapsed = time.time() - self.start_time
            sys.stderr.write("\r" + self.modname + self._percentage(prg) \
                             + " [" + barString + "] " + self._getTimeStr(prg, self.secondsElapsed))
            self.lastPercentage = self._percentage(prg)

    def _percentage(self, pbar):
        return '%3d%%' % (pbar * 100)

    def _formatTime(self, seconds):
        return str(int(seconds/3600)).zfill(2) + time.strftime(':%M:%S', time.gmtime(seconds))

    def _getTimeStr(self, prg=0, secondsElapsed=0):
        if prg == 0:
            return 'Time remaining: --:--:--'
        elif prg == 1:
            return 'Elapsed time:   %s\n' % self._formatTime(secondsElapsed)
        else:
            return 'Time remaining: %s' % self._formatTime(secondsElapsed / prg - secondsElapsed)

