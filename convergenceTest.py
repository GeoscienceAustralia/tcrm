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


import os
import io
import sys

import matplotlib
#matplotlib.use('Agg', warn=False)  # Use matplotlib backend

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
sns.set_context("paper")
figsize=(6.5, 4.5)
sns.set_style("whitegrid")


# Load the configuration file from the TCHA18, then open the database
# and get teh list of available locations.

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

def calculateARI(data, years):
    emprp = empReturnPeriod(np.sort(data))
    return np.sort(data)[-years:], emprp[-years:]

def bootstrap(data, n=1000, q=[5, 95], years=10000):
    d = np.empty((years, n))
    r = np.empty((years, n))
    for i in range(n):
        subset = np.random.choice(data, int(len(data/2)))
        d[:, i], r[:, i] = calculateARI(subset, years=10000)
    return np.percentile(d, q, axis=1), np.percentile(r, q, axis=1)

def plotConvergenceTest(locName):
    locId = locations['locId'][locations['locName']==locName][0]
    locLon = locations['locLon'][locations['locId']==locId][0]
    locLat = locations['locLat'][locations['locId']==locId][0]

    records = database.queries.locationRecords(db, str(locId))
    recs = records['wspd'][records['wspd'] > 0]
    data = np.zeros(int(NumSimulations*365.25))
    data[-len(recs):] = recs
    sortedmax = np.sort(data)
    emprp = empReturnPeriod(data)
    dd, rr = bootstrap(data, n=100)

    ep = 1./emprp
    ep1 = 1./rr[0,:]
    ep2 = 1./rr[1,:]
    mn = dd.mean(axis=0)
    delta = np.abs(np.diff(dd, axis=0))
    fdelta = delta/mn

    fig, ax1 = plt.subplots(1, 1, figsize=figsize)

    ax1.fill_between(rr[0,:], dd[1,:], dd[0,:], alpha=0.5, label="95th percentile")
    ax1.plot(emprp[-10000:], data[-10000:], color='k', label="Mean ARI")
    ax1.set_xscale('log')

    xlabel = 'Average recurrence interval (years)'
    ylabel = 'Wind speed (m/s)'
    title = "ARI wind speeds at " + locName + \
        " \n(%5.2f,%5.2f, n=%d)"%(locLon, locLat, len(recs))
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    ax1.legend(loc=2)
    addARIGrid(ax1)
    fig.tight_layout()
    plt.savefig(os.path.join(plotPath, "{0:05d}_ARI.png".format(locId)), 
                bbox_inches='tight')
    plt.close()

    fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
    ax2.semilogy(sortedmax[-10000:], ep[-10000:], color="k", label="Mean AEP")
    ax2.fill_betweenx(1./rr[0,:], dd[0,:], dd[1,:], alpha=0.5, label="95th percentile")
    ax2.set_xlabel(ylabel)
    title = "AEP wind speeds at " + locName + \
        " \n(%5.2f,%5.2f, n=%d)"%(locLon, locLat, len(recs))
    ax2.set_ylabel("Exceedance probability (events/year)")

    ax2.set_title(title)
    ax2.legend(loc=1)
    addAEPGrid(ax2)
    fig.tight_layout()
    plt.savefig(os.path.join(plotPath, "{0:05d}_AEP.png".format(locId)), 
                bbox_inches='tight')
    plt.close()
"""
    fig3, (ax3, ax4) = plt.subplots(2, 1, sharex=True)
    ax3.fill_between(emprp[emprp > 1][0:-1:2], fdelta,
                     color="#006983", alpha=0.5)
    ax3.set_ylabel('Fractional difference')
    ax3.set_title("Difference in convergence test ARI wind speeds at " + \
                  locName + " \n(%5.2f,%5.2f, n=%d)"%(locLon, locLat, len(recs)))
    ax3.set_xscale('log')
    addARIGrid(ax3)

    ax4.fill_between(emprp[emprp > 1][0:-1:2], delta,
                     color="#006983", alpha=0.5)
    ax4.set_ylabel("Difference (m/s)")
    ax4.set_xlabel(xlabel)
    ax4.set_xscale('log')
    addARIGrid(ax4)

    fig.tight_layout()
    plt.savefig(os.path.join(plotPath, "{0:05d}_ARI_delta.png".format(locId)), 
                bbox_inches='tight')
    plt.close()
"""
# Run the next cell, then select a location from the dropdown list and
# click the `"Run plotConvergenceTest"` button. This will take a
# minute or so to run, as it needs to extract all the values from the
# database, which is a significant number of events when there's
# 10,000 years of events.





locList = ['Carnarvon Airport', 'Port Hedland Airport',
           'Broome Airport', 'Darwin Airport',
           'Cairns Airport', 'Townsville Amo',
           'Rockhampton Airport', 'Willis Island']

#for locName in locList:
#    print(locName)
#   plotConvergenceTest(locName)

def plotConvergence(ax, locName):
    locId = locations['locId'][locations['locName']==locName][0]
    locLon = locations['locLon'][locations['locId']==locId][0]
    locLat = locations['locLat'][locations['locId']==locId][0]

    records = database.queries.locationRecords(db, str(locId))
    recs = records['wspd'][records['wspd'] > 0]
    data = np.zeros(int(NumSimulations*365.25))
    data[-len(recs):] = recs
    sortedmax = np.sort(data)
    emprp = empReturnPeriod(data)
    dd, rr = bootstrap(data, n=100)
    ax.plot(emprp[-10000:], data[-10000:], color='k', label="Mean ARI")
    ax.fill_between(rr[0,:], dd[1,:], dd[0,:], alpha=0.5, label="95th percentile")
    ax.set_xscale('log')
    #xlabel = 'Average recurrence interval (years)'
    #ylabel = 'Wind speed (m/s)'
    title = "{0} (n={1:d})".format(locName, len(recs))

    ax.set_title(title)
    addARIGrid(ax)


fig, axes = plt.subplots(4, 2, figsize=(6,8), sharex=True, sharey=True)
axlist = axes.flatten()
for i, loc in enumerate(locList):
    plotConvergence(axlist[i], loc)

axlist[0].legend(loc=2)
axlist[0].set_ylabel('Wind speed (m/s)')
axlist[2].set_ylabel('Wind speed (m/s)')
axlist[4].set_ylabel('Wind speed (m/s)')
axlist[6].set_ylabel('Wind speed (m/s)')

axlist[6].set_xlabel('Average recurrence interval (years)')
axlist[7].set_xlabel('Average recurrence interval (years)')

fig.tight_layout()
plt.savefig(os.path.join(plotPath, "ARI_convergence.png"),
                bbox_inches='tight')
