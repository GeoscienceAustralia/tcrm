.. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

.. _modelsetup:

====================
Setting up the model
====================

Execution of TCRM is controlled by reading the simulation settings
from a configuration file. The configuration file is a text file, and
can be edited in any text editor (e.g. Notepad, Wordpad, vi, emacs,
gedit). An example configuration file is provided in the examples
folder to give users a starting point.


.. _configurationfile:

The configuration file
======================

The TCRM configuration file is divided into a series of sections, each
with a set of option/value pairs. Most options have default values and
may not need to be specified in the configuration file. One value that
has no default is the Region gridLimit option. This defines the model
domain and must be set in any configuration file used.

.. _configureactions:

Actions
-------

This section defines which components of TCRM will be
executed. The options are:

* `DownloadData` - download input datasets (defaults are included)
* `DataProcess` - process the input TC track database
* `ExecuteStat` - calculate the TC statistics over the model domain
* `ExecuteTrackGenerator` - generate a set of stochastic TC tracks
* `ExecuteWindfield` - Calculate the wind field around a set of TC
  tracks
* `ExecuteHazard` - Calculate the return period wind speeds from a set
  of wind field files
* `PlotHazard` - Plot the return period wind speed maps and return
  period curves for locations in the model domain
* `PlotData` - Plot some basic statistical analyses of the input TC
  track database
* `ExecuteEvaluate` - Evaluate a set of stochastic TC tracks, comparing
  to the input TC track database.

All options are boolean (i.e. ``True`` or ``False``). ::

    [Actions]
    DataProcess = True
    ExecuteStat = True
    ExecuteTrackGenerator = True
    ExecuteWindfield = True
    ExecuteHazard = True
    PlotHazard = True
    PlotData = False
    ExecuteEvaluate = False
    CreateDatabase = True
    DownloadData = True

.. _configureregion:

Region
------

This section defines the simulation domain and the size of the grid over
which statistics are calculated. The simulation domain (``gridLimit``) is
specified as a Python dict with keys of ``xMin``, ``xMax``, ``yMin``
and ``yMax``. This sets the domain over which the wind fields and
hazard will be calculated. Stochastic tracks are generated over a
broader domain (called the "track domain"). The ``gridSpace`` option controls
the size of the grid cells, which are used for calculating statistics. At this
time, the values here must be integer values, but can be different in the ``x``
(east-west) and ``y`` (north-south) directions. The ``gridInc`` option
control the incremental increase in grid cell size when insufficient
observations are located within a grid cell (see the :mod:`StatInterface`
description)::

    [Region]
    gridLimit = {'xMin': 113.0, 'xMax': 124.0, 'yMin': -24.0, 'yMax': -13.0}
    gridSpace = {'x':1.0,'y':1.0} 
    gridInc = {'x':1.0,'y':0.5}

.. _configuredataprocess:

DataProcess
-----------

This section controls aspects of the processing of the input track
database. Firstly, the ``InputFile`` option specifies the file to be
processed. A relative or absolute path can be used. If no path name is
included (as in the example below), then TCRM assumes the file is
stored in the ``input`` path. If using an automatically
downloaded dataset, then this file name must match the name
specified in the appropriate dataset section (which is named by the
``Source`` option in this section) of the configuration file (further
details below).

The ``Source`` option is a string value that acts as a pointer to a
subsequent section in the configuration file, that holds details of
the input track file structure. The additional section must have the same label
as set here (the label is case-sensitive).

The ``StartSeason`` and ``FilterSeason`` options control what years of
the input track database are used in calibrating the model. In the
default case, only data from 1981 onwards is used for model
calibration. If ``FilterSeasons = False``, no season filtering is
performed and the full input track database is used. ::

    [DataProcess]
    InputFile = ibtracs.since1980.list.v04r00.csv
    StartSeason = 1981
    FilterSeasons = True
    Source = IBTRACS

.. _configurestatinterface:

StatInterface
-------------

The ``StatInterface`` section controls the methods used to calculate
distributions of TC parameters from the input track database.

``kdeType`` specifies the kernel used in the kernel
density estimation method for creating probability density functions
that are used in selecting initial values for the stochastic TC events
(e.g. longitude, latitude, initial pressure, speed and
bearing). ``kdeStep`` defines the increment in the generated
probability density functions and cumulative distribution functions.

Options for ``kdeType`` ::
    'gau'
    'epa'
    'uni'
    'tri'
    'biw'
    'triw'
    'cos'
    'cos2'


``kde2DType`` is deprecated.

``minSamplesCell`` sets the minimum number of valid observations in
each grid cell that are required for calculating the distributions,
variances and autocorrelations used in the :mod:`TrackGenerator`
module. If there are insufficient valid observations, then the bounds
of the grid cell are incrementally increased (in steps as specified by
the ``gridInc`` values) until sufficient observations are found. ::

    [StatInterface]
    kdeType = gau
    kde2DType = gau
    kdeStep = 0.2
    minSamplesCell = 100

.. _configuretrackgenerator:

TrackGenerator
--------------

The ``TrackGenerator`` section controls the stochastic track
generation module. It is here that users can control the number of
events and the number of years generated.

The ``NumSimulations`` option sets the number of TC event sets that
will be generated. Any integer number of events (up to 1,000,000) is
possible. ``YearsPerSimulation`` sets the number of simulated years
that will be generated for each event set. For evaluating hazard, the
value should be set to 1, as the extreme value distribution fitting
process assumes annual maxima. The annual frequency of events is based
on a Poisson distribution around the mean annual frequency, which is
determined from the input track database.

For track model evaluations, it is recommended to set
``YearsPerSimulation`` to a similar number to the number of years in
the input track database. For example, in our testing that used data
from 1981--2013, we set the value to 30.

``NumTimeSteps`` controls the maximum lifetime an event can exist
for. ``TimeStep`` sets the time interval (in hours) for the track
generator. 

``SeasonSeed`` and ``TrackSeed`` are used to fix the random
number generator on parallel systems to ensure truly random numbers on
each individual processor. If they are absent, the seed is set using an integer
representation of the current time, and is recorded in the output metadata (e.g.
attributes in the netcdf files). ::

    [TrackGenerator]
    NumSimulations = 500
    YearsPerSimulation = 1
    NumTimeSteps = 360
    TimeStep = 1.0
    SeasonSeed = 1
    TrackSeed = 1

This example will generate 500 realisations of one year of TC activity, with
hourly timesteps to a maximum of 360 hours. 


.. _configurewindfield:

WindfieldInterface
------------------

The ``WindfieldInterface`` section controls how the wind fields from
each track in the simulated tracks are calculated. There are two main
components to the wind field -- the radial profile and the boundary
layer model.

The ``profileType`` option sets the radial profile used. Valid values are:

* ``holland`` -- the radial profile of Holland (1980) [1]_
* ``powell`` -- Similar to the Holland profile, but uses a variable
  beta parameter that is a function of latitude and size. [2]_
* ``schloemer`` -- From Schloemer (1954) -- essentially the Holland
  profile with a beta value of 1 [3]_
* ``willoughby`` -- From Willoughby and Rahn (2004). Again, the
  Holland profile, with beta a function of the maximum wind speed,
  radius to maximum wind and latitude [4]_
* ``jelesnianski`` -- From Jelesnianski (1966). [5]_
* ``doubleHolland`` -- A double exponential profile from McConochie
  *et al.* (2004) [6]_

The ``windFieldType`` value selects the boundary layer model
used. Three boundary layer models have been implemented:

* ``kepert`` -- the linearised boundary layer model of Kepert (2001)
  [7]_
* ``hubbert`` -- a vector addition of forward speed and tangential
  wind speed from Hubbert *et al.* (1994) [8]_
* ``mcconochie`` -- a second vector addition model, from McConochie
  *et al.* (2004) [6]_

The ``beta`` option specifies the |beta| parameter used in the Holland
wind profile. The additional |beta| options (``beta1`` and ``beta2``)
are used in the ``doubleHolland`` wind profile, which is a double
exponential profile, therefore requiring two |beta| parameters.

``thetaMax`` is used in the McConochie and Hubbert boundary layer
models to specify the azimuthal location of the maximum wind speed
under the translating storm.

``Margin`` defines the spatial extent over which the wind field is
calculated and is in units of degrees. A margin of 5 is recommended
for hazard models, to ensure low wind speeds from distant TCs are
incorporated into the fitting procedure.

``Resolution`` is the horizontal resolution (in degrees) of the wind
fields. Values should be no larger than 0.05 degrees, as the absolute
peak of the radial profile may not be adequately resolved, leading to
an underestimation of the maximum wind speeds. 

``Domain`` is an alternative to setting ``Margin``. If set to "bounded"
(default), the wind field domain will be determined by the ``Margin`` option. 
If set to "full", the wind field domain will be set to match ``Region -
gridLimit``. 

*If this option is chosen, the execution time will be significantly
longer. Recommended for single scenarios only*. ::

    [WindfieldInterface]
    profileType = holland
    windFieldType = kepert
    beta = 1.3
    beta1 = 1.3
    beta2 = 1.3
    thetaMax = 70.0
    Margin = 2
    Resolution = 0.05
    Domain = bounded

.. _configurehazard:

Hazard
------

The ``Hazard`` section controls how the model calculates average recurrence
(ARI) wind speeds, and whether to calculate confidence ranges.

``ExtremeValueDistribution`` sets the method for calculating ARI wind speeds.
Options are "emp" (empirical), "power" (power law), "GPD" (Generalised pareto
distribution) or "GEV" (Generalised Extreme Value distribution). 

The ``Years`` option is a comma separated list of integer values that
specifies the return periods for which wind speeds will be
calculated. For ``ExtremeValueDistribution = emp``, the years cannot exceed the
total number of simulated years in the ``TrackGenerator`` options. 

``MinimumRecords`` sets the minimum number of values
required for performing the fitting procedure at a given grid point.

``CalculateCI`` sets whether the :mod:`hazard` module will calculate
confidence ranges using a bootstrap resampling method. If ``True``,
the module will run the fitting process multiple times and calculate
upper and lower percentile values of the resulting return period wind
speeds. The ``PercentileRange`` option sets the range -- for a value
of 90, the module will calculatae the 5th and 95th percentile
values. ``SampleSize`` sets the number of randomly selected values
that will be used in each realisation of the extreme value fitting
procedure for calculating the confidence range. 

``SmoothPlots`` will apply a gaussian filter to the data before plotting on maps
to minimise the inference of lines on maps. This may cause the maps to have
large areas of no data due to the filtering function. ::

    [Hazard]
    ExtremeValueDistribution = emp
    Years = 2,5,10,20,25,50,100,200,250,500,1000
    MinimumRecords = 50
    CalculateCI = True
    PercentileRange = 90
    SampleSize = 50
    PlotSpeedUnits = mps
    SmoothPlots = True

.. _configurermw:

RMW
----

The ``RMW`` section contains a single option: ``GetRMWDistFromInputData``. 
Set this value to ``True`` if the input track database has reliable data 
on the radius to maximum winds. 

If no suitable data exists (``GetRMWDistFromInputData = False``), TCRM will use
a regression model to determine RMW from the intensity and latitude of the
storm. ::

    [RMW]
    GetRMWDistFromInputData = False

.. _configureinput:

Input
-----

The ``Input`` section sets the source of some supplementary data, as
well as the datasets to be automatically downloaded. The ``LandMask``
option specifies the path to a netcdf file (supplied) that contains a
land/sea mask. The ``MSLPFile`` option specifies the path to a netcdf
file (downloaded) that contains daily long-term mean sea level
pressure data (e.g. from a NCEP/NCAR reanalysis products). The 
``LocationFile`` option specifies the path to a point shape file that contains
the longitude and latitude of locations for which to extract hazard information
at the conclusion of a simulation.

The ``Datasest`` option is a comma separated list of values indicating
the data that should be downloaded on first execution. For each value
in the list, there must be a corresponding section in the
configuration file, that has options of ``URL`` (the URL of the data
to be downloaded), ``path`` (where to store the data once it has been
downloaded) and ``filename`` (the filename to give to the data once
downloaded).

In the example below, for the ``IBTRACS`` dataset, there are
additional options that describe the format of the track database with
the same name.  This is a legitimate approach, so long as there are no
duplicate options.

Note that the ``filename`` option in the ``IBTRACS`` section matches
the ``InputFile`` option in the ``DataProcess`` section, and the
``filename`` in the ``LTMSLP`` section matches the ``MSLPFile`` in the
``Input`` section.

The ``CoastlineGates`` option specifies the path to a comma-delimited
text file that holds the points of a series of coastline gates that
are used in the :mod:`Evaluate.landfallRates` module. ::

    [Input]
    LocationFile = input/stationlist.shp
    LandMask = input/landmask.nc
    MSLPFile = MSLP/slp.day.ltm.nc
    Datasets = IBTRACS,LTMSLP
    CoastlineGates = input/gates.csv

    [IBTRACS]
    URL = https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.since1980.list.v04r00.csv
    path = input
    Filename = ibtracs.since1980.list.v04r00.csv
    Columns = tcserialno,season,num,skip,skip,skip,date,skip,lat,lon,skip,pressure
    FieldDelimiter = ,
    NumberOfHeadingLines = 2
    PressureUnits = hPa
    LengthUnits = km
    SpeedUnits = kph
    DateFormat = %Y-%m-%d %H:%M:%S

    [LTMSLP]
    URL = ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/slp.day.1981-2010.ltm.nc
    path = MSLP
    filename = slp.day.ltm.nc

.. _configureoutput:

Output
------

The ``Output`` section defines the destination of the model output. Set the 
``Path`` option to the directory where you wish to store the data. Paths can 
be relative or absolute. By default, output is stored in a subdirectory of 
the working directory named ``output``. ::

    [Output]
    Path = output

.. _configurelogging:

Logging
-------

The ``Logging`` section controls how the model records progress to
file (and optionally STDOUT). ``LogFile`` option specifies the name of
the log file. If no path is given, then the log file will be stored in
the current working directory. For parallel execution, a separate log
file is created for each thread, with the rank of the process appended
to the name of the file.

The ``LogLevel`` is one of the :mod:`Logging` `levels
<https://docs.python.org/2/library/logging.html#logging-levels>`_. Default
is ``INFO``. 

The ``Verbose`` option allows users to print all logging
messages to the standard output (default False). This can be useful when attempting to
identify problems with execution. For parallel execution, this is set
to ``False`` (to prevent repeated messages being printed to the
screen). 

Setting the ``ProgressBar`` option to ``True`` will display a
simple progress bar on the screen to indicate the status of the model
execution (default False). This will be turned off if TCRM is executed on a parallel
system, or if it is run in batch mode. 

If ``Datestamp = True``, a timestamp will be included in the filename for the
log file (default False). ::

    [Logging]
    LogFile = main.log
    LogLevel = INFO
    Verbose = False
    ProgressBar = False
    Datestamp = False

.. _configuresource:

Source format options
---------------------

For the input data source specified in the :menuselection:`DataProcess --> Source`
option, there must be a corresponding section of the given name. In
this example case, the source is specified as ``IBTRACS`` (the same as
one of the ``Dataset`` options). The ``IBTRACS`` section therefore
controls both the download dataset options, and specifies the textural
format of the input track database.

The options that relate to the dataset download are ``URL``, ``path``
and ``filename``. ``URL`` specifies the location of the data to be
downloaded. The ``path`` option specifies the path name for the
storage location of the dataset. The ``filename`` option gives the
name of the file to be saved (this can be different from the name of
the dataset).

The remaining options relate to the format of the track
database. ``Columns`` is a comma-separated list of the column names in
the input database. If a column is to be ignored, it should be named
``skip``. The ``FieldDelimiter`` is the delimiter used in the input
track database (it's assumed that the input file is a text format
file!). The ``NumberOfHeadingLines`` indicates the number of text
lines at the top of the file that should be ignored (usually this is
column headers -- due to the multiple lines used in some track
databases, TCRM does not attempt to decipher the column names from the
header. ``PressureUnits``, ``LengthUnits`` and ``SpeedUnits`` specify
the units the numerical values of pressure, distance and speed
(respectively) used in the input track database. The ``DateFormat``
option is a string represenation of the date format used in the track
database. The format should use Python's `datetime
<https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior>`_
formats.  ::

    [IBTRACS]
    URL = https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.since1980.list.v04r00.csv
    path = input
    Filename = ibtracs.since1980.list.v04r00.csv
    Columns = tcserialno,season,num,skip,skip,skip,date,skip,lat,lon,skip,pressure
    FieldDelimiter = ,
    NumberOfHeadingLines = 2
    PressureUnits = hPa
    LengthUnits = km
    SpeedUnits = kph
    DateFormat = %Y-%m-%d %H:%M:%S
 
.. _references:

References
----------

.. [1] Holland, G. J. (1980): An Analytic Model of the Wind and Pressure 
       Profiles in Hurricanes. *Monthly Weather Review*, **108**
.. [2] Powell, M., G. Soukup, S. Cocke, S. Gulati, N. Morisseau-Leroy, S. 
       Hamid, N. Dorst, and L. Axe (2005): State of Florida hurricane loss 
       projection model: Atmospheric science component. *Journal of Wind 
       Engineering and Industrial Aerodynamics*, **93**, 651--674
.. [3] Schloemer, R. W. (1954): Analysis and synthesis of hurricane wind 
       patterns over Lake Okeechobee. *NOAA Hydrometeorology Report* **31**, 
       1954
.. [4] Willoughby, H. E. and M. E. Rahn (2004): Parametric Representation 
       of the Primary Hurricane Vortex. Part I: Observations and 
       Evaluation of the Holland (1980) Model. *Monthly Weather Review*, 
       **132**, 3033--3048
.. [5] Jelesnianski, C. P. (1966): Numerical Computations of Storm Surges 
       without Bottom Stress. *Monthly Weather Review*, **94**, 379--394
.. [6] McConochie, J. D., T. A. Hardy, and L. B.  Mason (2004):  Modelling 
       tropical cyclone over-water wind and pressure fields. *Ocean 
       Engineering*, **31**, 1757--1782

.. [7] Kepert, J. D. (2001): The Dynamics of Boundary Layer Jets 
       within the Tropical Cyclone Core. Part I: Linear Theory.  
       *J. Atmos. Sci.*, **58**, 2469--2484 
.. [8] Hubbert, G. D., G. J. Holland, L. M. Leslie and M. J. Manton (1991): 
       A Real-Time System for Forecasting Tropical Cyclone Storm Surges. 
       *Weather and Forecasting*, **6**, 86--97

