.. _scenariomodelling:
Scenario modelling
==================

TCRM can also be used to simulate the wind field from an individual
event. Given the track of a tropical cyclone, users can run the
:mod:`wind` module only, and generate the maximum wind swath from a
TC. The `tcevent.py` script enables users to efficiently run a
scenario simulation, including performing temporal interpolation of
the track positions to generate a realistic representation of the wind
field. This can be useful for exploring synthetic events in more
detail, or for analysing the wind field from an historical tropical
cyclone.

TCRM does not currently account for local landscape effects on wind
speed (e.g. topographic enhancement or changes in surface roughness
due to vegetation or built environments). Users should determine the
best way to include these effects for their own purposes.

As TCRM uses a parametric profile, the primary vortex is axisymmetric,
and asymmetry in the surface winds arises solely due to the forward
motion of the cyclone and the (uniform) surface friction. For a
complete and accurate representation of the winds from an actual
tropical cyclone, the simulated winds from TCRM should be combined
with observed maximum wind speeds in the vicinity of the cyclone using
spatial interpolation methods (e.g. kriging).

.. _scenariosetup:
Setting up a scenario
---------------------

Users will need to manually edit a configuration file to execute a
scenario simulation. Many of the options in the main TCRM
configuration file are not required for scenario simulation, so it is
recommended to create a reduced configuration file that contains only
the required sections and options.

A basic configuration file for a scenario simulation would look like
this::

    [DataProcess]
    InputFile=scenario.csv

    [WindfieldInterface]
    Source=NRL    # Naval Research Laboratory
    Margin=3
    Resolution=0.02
    profileType=powell
    beta=1.6
    windFieldType=kepert
    TimeStep=0.0833333333

    [Timeseries]
    Extract=True
    StationFile=./input/stationlist.shp
    StationID=WMO

    [Output]
    Path = ./output/scenario
    
    [NRL]
    Columns=index,skip,skip,date,lat,lon,skip,skip,pressure,rmax
    FieldDelimiter=,
    NumberOfHeadingLines=0
    PressureUnits=hPa
    SpeedUnits=kts
    LengthUnits=km
    DateFormat=%Y%m%d %H:%M



.. _runningscenario:
Running the scenario
--------------------

The `tcevent.py` script loads the track file, performs the temporal
interpolation and then passes the interpolated track file to the
:mod:`wind` module. Since the ``Region`` section has been removed, the
:meth:`wind.run` method will set the grid domain to cover the entire
extent of the track. 

Make sure ``python`` is in your system path, then from the base
directory, call the ``tcevent.py`` script, with the configuration 
file option included. For example, to run the example scenario::

    python tcevent.py -c example/scenario.ini

If the ``-v`` option is included, then all logging messages will be
printed to the console. The level of logging detail is set in the
configuration file.

:note: `tcevent.py` cannot be executed in parallel. 

.. _scenariocmdlineargs:
Command line arguments
----------------------

 -c file, --config file   Path to a configuration file.
 -v, --verbose            If given, logging messages will be printed 
                          to the console.
 -d, --debug              In the case that execution results in an exception, 
                          allow the Python stack to call into the stack trace 
                          (through implementation of a custom hook script) and 
                          start the Python debugger (:mod:`pdb`). 

.. _timeseries:
Extract time series data
------------------------

When running a scenario, it is possible to extract a time series of
the wind speed and sea level pressure values from the grid at selected
locations. The locations are defined in a user-supplied point shape
file (a location database is planned for inclusion in future
versions to better facilitate this feature). The shape file should
contain a field with a unique identifier, otherwise station output is
numbered sequentially through the locations.

* Locations must be provided in geographic coordinates (longitude,
  latitude coordinates). No reprojection is performed.
* Output is a 'regional' wind speed -- that is, the wind speed at that
  location, excluding local topographic or landscape effects. These
  effects can be incorporated offline (i.e. outside the TCRM
  framework).

The data is stored in a separate csv file for each location, and data
is plotted on a simple figure for visual inspection.

.. figure:: ./timeseries_example.png
    :align: center
    :alt: Time series example from Guiuan, Philippines
    :figclass: align-center

    Time series data for Super Typhoon *Haiyan* at Guiuan, Samar, Philippines. 

:Note: The double labels on the secondary (right-hand) y-axis require
       Matplotlib version 1.3 or later.

