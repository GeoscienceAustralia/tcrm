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


Title: trackSimulation.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2007-07-27
Description: Generate a set of events and generate landfall statistics
             for each set.
Reference:

SeeAlso:

Constraints:

Version: 106
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-01-13
Modification: Removed landfall PDF calculations

Version: $Rev: 643 $
ModifiedBy:
ModifiedDate:
Modification:

$Id: trackSimulation.py 643 2011-10-31 05:32:34Z nsummons $
"""


import os, sys, pdb, logging

from numpy import random
import math
import TrackGenerator as TG
from Utilities.progressbar import ProgressBar

from Utilities.config import cnfGetIniValue
#from Utilities.files import flSaveFile

__version__ = '$Id: trackSimulation.py 643 2011-10-31 05:32:34Z nsummons $'

def trackSimulation(configFile, nSim, yrsPerSim, meanFrequency, outputPath, format='csv',
                    dt=1, tsteps=360, maxDist=2500, autoCalc_gridLimit=None):
    """
    Simulate nSim sets of cyclones, with a mean frequency and a given
    number of years in each simulation
    """
    logger = logging.getLogger()
    pbar = ProgressBar('(3/6) Generating tracks:     ')
    numCyclones = sum(random.poisson(meanFrequency, yrsPerSim))
    tracks = TG.TrackGenerator(configFile, dt, tsteps, autoCalc_gridLimit=autoCalc_gridLimit, progressbar=pbar)

    for i in range(1,nSim+1):
        # This seemingly simple expression provides a method to generate
        # a timeseries of annual TC numbers with a given mean.
        numCyclones = sum(random.poisson(meanFrequency, yrsPerSim))
        logger.info("Generating %d synthetic TC events for simulation %04d"
                     % (numCyclones,i))
        trackFile = os.path.join(outputPath,"tracks.%04d.%s"%(i,format))
        tracks.generatePath(nEvents=numCyclones, outputFile=trackFile)
        pbar.update(i/float(nSim), 0.3, 1.0)
        #flSaveFile(trackFile, trackData,  delimiter=',', fmt='%7.2f')
        #data = transpose(array([trackData[:,0], trackData[:,2], trackData[:,3], trackData[:,6]]))

        #LF.loadData(data)
        #LF.crossing(outputFile=countFile)
        #pdfInd,pdfLand,pdfOff = LF.generatePDF()
        #if i == 1:
        #    pdfLandfall = [transpose(array(pdfInd)), transpose(array(pdfLand))]
        #    pdfOffshore = [transpose(array(pdfInd)), transpose(array(pdfOff))]
        #else:
        #    pdfLandfall = concatenate((pdfLandfall,array([pdfLand])))
        #    pdfOffshore = concatenate((pdfOffshore,array([pdfOff])))

        #print "Completed set %i of %04d event sets"%(i,numSimulations)

    #flSaveFile(os.path.join(outputPath,"summary_landfallpdf.dat"),transpose(pdfLandfall),header='',delimiter=',',fmt='%9.5f')

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg
        
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s',
                        filename=cnfGetIniValue(configFile, 'Logging', 'LogFile', __file__.rstrip('.py') + '.log'),
                        filemode='w')
    numSimulations = cnfGetIniValue(configFile, 'TrackGenerator',
                                    'NumSimulations', 50)
    dt = cnfGetIniValue(configFile, 'TrackGenerator', 'TimeStep', 1)
    tsteps = cnfGetIniValue(configFile, 'TrackGenerator', 'NumTimeSteps', 360)

    gridLimit = eval(cnfGetIniValue(configFile, 'TrackGenerator',
                                    'gridLimit'))
    outputPath = os.path.join(cnfGetIniValue(configFile, 'Output', 'Path'),
                              'tracks')
    yrsPerSim = cnfGetIniValue(configFile, 'TrackGenerator', 'YearsPerSimulation', 10)
    frequency = cnfGetIniValue(configFile, 'TrackGenerator', 'Frequency')

    trackSimulation(configFile, numSimulations, frequency, yrsPerSim, gridLimit,
                    outputPath, dt=dt, tsteps=tsteps)

