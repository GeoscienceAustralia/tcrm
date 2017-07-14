Evaluation
==========

TCRM provides functionality to evaluate the performance of the track
model over the model domain through the :mod:`Evaluate` module. By
running a number of simulations, each generating a similar number of
years of events to the input track database, we can compare the
synthetic events with the observations. Metrics include track density,
genesis distribution, minimum pressure distributions, landfall rates
and longitude crossing rates. 


Setup
-----

An example configuration for executing an evaluation. Other sections
of the configuration file should remain unchanged. ::

    [Actions]
    DataProcess=True
    ExecuteStat=True
    ExecuteTrackGenerator=True
    ExecuteWindfield=False
    ExecuteHazard=False
    PlotHazard=False
    PlotData=False
    ExecuteEvaluate=True
    DownloadData=True

    [Input]
    LandMask=input/landmask.nc
    MSLPFile=MSLP/slp.day.ltm.nc
    CoastlineGates=input/gates.txt

    [TrackGenerator]
    NumSimulations=500
    YearsPerSimulation=30
    NumTimeSteps=360
    TimeStep=1.0
    SeasonSeed=1
    TrackSeed=1

Note that the ``YearsPerSimulation`` option is now set to 30, so TCRM
will generate 500 event sets, each with the equivalent of 30 years of
TCs. The annual frequency of events is based on a Poisson distribution
around the mean annual frequency, which is determined from the input
track database.

A csv-format file containing a set of gates around the coastline is
required for :mod:`Evaluate.landfallRates`. Gates are defined as a
series of points around the coastline in number, longitude, latitude
format. Below are the first 15 points in the sample gates provided in
the code base. ::

      0.00, 114.78, -34.70
      1.00, 115.20, -32.95
      2.00, 114.90, -31.17
      3.00, 114.49, -29.42
      4.00, 113.66, -27.83
      5.00, 112.73, -26.28
      6.00, 112.96, -24.50
      7.00, 113.22, -22.72
      8.00, 114.47, -21.43
      9.00, 115.98, -20.44
      10.00, 117.74, -20.08
      11.00, 119.46, -19.55
      12.00, 121.12, -18.84
      13.00, 121.74, -17.15
      14.00, 123.04, -15.91
      15.00, 124.49, -14.85

Output
------

The :mod:`Evaluate` module generates a set of figures that compare the
synthetic event set to the input track database. For the synthetic
events, the module calculates a mean value and upper and lower
percentiles of the distribution (commonly the 5th and 95th percentile
values). 

Some of the results are also saved to netCDF files (track density,
pressure distributions and longitude crossing rates) for users to
import into other graphics packages.

