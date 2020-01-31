"""
:mod:`progressbar` -- display progress on STDOUT
================================================

.. module:: progressbar
    :synopsis: print a progress bar on to standard output.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import time
import sys
import logging

log = logging.getLogger()

class ProgressBar(object):

    def __init__(self, modname, showbar=True):
        self.modname = modname #+ " "
        self.showbar = False
        self.lastPercentage = None
        self.screenWidth = 79
        self.barWidth = (self.screenWidth - len(self.modname) -
                         len(self._getTimeStr()) - 8)
        self.start_time = time.time()
        self.secondsElapsed = 0
        self.update(0.0)
        if (sys.stderr.isatty() and sys.stdin.isatty()):
            self.showbar = showbar

    def update(self, progress, startPos=0, endPos=1):
        if self.showbar:
            prg = progress * (endPos - startPos) + startPos
            if self._percentage(prg) != self.lastPercentage:
                barfill = int(round(self.barWidth * prg))
                barString = (''.join(['#' for i in range(barfill)]
                             + [' ' for i in range(self.barWidth - barfill)]))
                self.secondsElapsed = time.time() - self.start_time
                message = "\r{0} {1} [{2}] {3}".\
                          format(self.modname,
                                 self._percentage(prg),
                                 barString,
                                 self._getTimeStr(prg, self.secondsElapsed))
                sys.stderr.write(message)

                #sys.stderr.write("\r" + self.modname + self._percentage(prg) \
                #                 + " [" + barString + "] "  \
                #                 + self._getTimeStr(prg, self.secondsElapsed))
                self.lastPercentage = self._percentage(prg)

    def _percentage(self, pbar):
        return '%3d%%' % (pbar * 100)

    def _formatTime(self, seconds):
        return (str(int(seconds/3600)).zfill(2)
                + time.strftime(':%M:%S', time.gmtime(seconds)))

    def _getTimeStr(self, prg=0, secondsElapsed=0):
        if prg == 0:
            return 'Time remaining: --:--:--'
        elif prg == 1:
            return 'Elapsed time:   %s\n' % self._formatTime(secondsElapsed)
        else:
            return ('Time remaining: %s' % \
                    self._formatTime(secondsElapsed / prg - secondsElapsed))


class SimpleProgressBar(ProgressBar):

    def __init__(self, modname, showbar=True):
        ProgressBar.__init__(self, modname, showbar)
        self.lastPercentage = 0.

    def update(self, progress, startPos=0, endPos=1, incr=5.):
        prg = progress * (endPos - startPos) + startPos
        percent = prg * 100.
        if (self.showbar and sys.stderr.isatty()):
            # throttle output
            if percent >= 99. or percent >= (self.lastPercentage + incr):
                message = "{0} {1} complete".\
                          format(self.modname,
                                 self._percentage(prg).strip().rjust(4))

                print(('********* ' + message), file=sys.stderr)
                log.info(message)
                self.lastPercentage += incr
