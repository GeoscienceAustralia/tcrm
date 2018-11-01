# # Empirical recurrence interval uncertainty
# 
# The TCHA provides a catalogue of 10,000 years of TC activity. This
# can be considered a large number of observations, compared to the
# average recurrence intervals (ARIs) that are normally considered for
# TC-related winds (e.g. 100-, 200-, 500-year ARIs).
# 
# Therefore, we can use the simulated wind speed values to directly
# estimate the ARI wind speeds, without needing to fit extreme value
# distributions - something that is necessary when using observed wind
# speeds with record lengths of only around 50 years.
# 
# Here, we calculate the empirical ARI wind speeds from the synthetic
# catalogue, but do so in a manner that explores the range of values
# at these ARIs.
# 
# Essentially, we subsample the 10,000 year catalogue into smaller
# collections, then calculate the empirical ARIs from these smaller
# collections. The range of ARI wind speed values from these
# subsampled collections then gives us a guide on the convergence of
# these values in the full collection - if the ARIs show only a small
# range, then the ARI values have converged. We would expect the ARIs
# to diverge at close to the total number of simulated years
# (i.e. close to 10,000 years).
# 

from __future__ import print_function, division
import os
import io
import sys

import matplotlib
matplotlib.use('Agg', warn=False)  # Use matplotlib backend

import database
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.ticker import LogLocator, FormatStrFormatter

from Utilities.config import ConfigParser

from extremes import returnLevels, empReturnPeriod, returnPeriodUncertainty, gpdSelectThreshold
from distributions import fittedPDF

import random

import seaborn as sns
sns.set_context("notebook")
sns.set_style("whitegrid")


# Load the configuration file from the TCHA18, then open the database
# and get teh list of available locations.

# In[2]:

configFile = "/home/547/cxa547/tcrmconfig/tcrm2.1.ini"
config = ConfigParser()
config.read(configFile)
outputPath = config.get('Output', 'Path')
plotPath = os.path.join(outputPath, 'plots','convergence')
NumSimulations = config.getint('TrackGenerator', 'NumSimulations')

db = database.HazardDatabase(configFile)
locations = db.getLocations()
locNameList = list(locations['locName'])

# The following step performs the calculations. First a helper
# function to add nicely formatted grid lines on a logarithmic axis.
# 
# The second function (`plotConvergenceTest`) loads the data from the
# database, then splits into two separate collections (called `d1` and
# `d2`). For each of these, we then calculate empirical ARI values and
# plot alongside each other. We also plot the annual exceedance
# probability as an alternate view on the likelihood of extreme winds.


def addARIGrid(axes):
    """
    Add a logarithmic graticule to the subplot axes.
    :param axes: :class:`axes` instance.
    """

    axes.xaxis.set_major_locator(LogLocator())
    axes.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    axes.xaxis.set_minor_locator(LogLocator(subs=[.1, .2, .3, .4, .5, .6, .7, .8, .9]))
    axes.autoscale(True, axis='x', tight=True)
    axes.grid(True, which='major', linestyle='-')
    axes.grid(True, which='minor', linestyle='--', linewidth=0.5)
    
def addAEPGrid(axes):
    """
    Add a logarithmic graticuyle to the subplot axes
    :param axes: :class:`axes` instance
    """
    axes.yaxis.set_major_locator(LogLocator())
    axes.yaxis.set_minor_locator(LogLocator(subs=[.1, .2, .3, .4, .5, .6, .7, .8, .9]))
    axes.autoscale(True, axis='y', tight=True)
    axes.grid(True, which='major', linestyle='-')
    axes.grid(True, which='minor', linestyle='--', linewidth=0.5)
    
    
def plotConvergenceTest(locName):
    locId = locations['locId'][locations['locName']==locName][0]
    locLon = locations['locLon'][locations['locId']==locId][0]
    locLat = locations['locLat'][locations['locId']==locId][0]

    records = database.locationRecords(db, locId)
    recs = records['wspd'][records['wspd'] > 0]
    data = np.zeros(int(NumSimulations*365.25))
    data[-len(recs):] = recs
    sortedmax = np.sort(data)
    emprp = empReturnPeriod(data)
    random.shuffle(data)
    d1 = data[:int(len(data)/2)]
    d2 = data[int(len(data)/2+1):]
    sortedmax1 = np.sort(d1)
    sortedmax2 = np.sort(d2)
    emprp1 = empReturnPeriod(d1)
    emprp2 = empReturnPeriod(d2)
    ep = 1./emprp
    ep1 = 1./emprp1
    ep2 = 1./emprp2
    
    fig, ax1 = plt.subplots(1, 1)
    ax1.semilogx(emprp[emprp > 1], sortedmax[emprp > 1], color='k', 
                 label="Mean ARI")
    ax1.semilogx(emprp2[emprp2> 1], sortedmax2[emprp2 > 1], color="#006983",
                 label="Convergence check 1")
    ax1.semilogx(emprp1[emprp1> 1], sortedmax1[emprp1 > 1], color="#A33F1F",
                 label="Convergence check 2")
    ax1.set_xscale('log')

    xlabel = 'Average recurrence interval (years)'
    ylabel = 'Wind speed (m/s)'
    title = "ARI wind speeds at " + locName + \
        ", \n(%5.2f,%5.2f, n=%d)"%(locLon, locLat, len(recs))
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    addARIGrid(ax1)
    fig.tight_layout()
    plt.savefig(os.path.join(plotPath, "{0:05d}_ARI.png".format(locId)), 
                bbox_inches='tight')
    plt.close()
    fig2, ax2 = plt.subplots(1, 1)
    ax2.semilogy(sortedmax[emprp > 1], ep[emprp > 1], color="k",
                 label="Mean exceedance rate")

    ax2.semilogy(sortedmax1[emprp1 > 1], ep1[emprp1 > 1], color="#006983",
                 label="Convergence check 1")
    ax2.semilogy(sortedmax2[emprp2 > 1], ep2[emprp2 > 1], color="#A33F1F",
                 label="Convergence check 2")
    ax2.set_xlabel(ylabel)
    title = "AEP wind speeds at " + locName + \
        ", \n(%5.2f,%5.2f, n=%d)"%(locLon, locLat, len(recs))
    ax2.set_ylabel("Exceedance probability")

    ax2.set_title(title)
    addAEPGrid(ax2)
    fig.tight_layout()
    plt.savefig(os.path.join(plotPath, "{0:05d}_AEP.png".format(locId)), 
                bbox_inches='tight')
    plt.close()

# Run the next cell, then select a location from the dropdown list and
# click the `"Run plotConvergenceTest"` button. This will take a
# minute or so to run, as it needs to extract all the values from the
# database, which is a significant number of events when there's
# 10,000 years of events.


for locName in locNameList:
    #print(locName)
    plotConvergenceTest(locName)
