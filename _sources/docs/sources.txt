
.. _sourceformats:

==============================================
Using different TC databases - source formats 
==============================================

There are a wide variety of formats of tropical cyclone track
databases available, both observational and model generated (e.g. we
have used tropical cyclone-like vortices derived from regional and
global climate models). To introduce users to the model, we include
the format for the `International Best Tracks Archive for Climate
Stewardship <http://www.ncdc.noaa.gov/oa/ibtracs/>`_ (IBTrACS-WMO
format). The default action for TCRM is to download the current
version (v03r05 at August 2014) and use this to calibrate the model.

.. _requiredfields:

Required fields
---------------
As a bare minimum, TCRM requires data on the date/time, longitude,
latitude, central pressure of each cyclone observation. It also
requires a field that indicates each unique TC event in the
database. This can be a ``tcserialno`` - a unique string for each
event; a TC ``num`` which is an integer value for each event in a
season; or an ``indicator`` field, which is set to ``1`` for the first
observation of a new TC, and ``0`` otherwise.

For individual scenario modelling, an additional field containing the
radius to maximum winds (``rmax``) must be included.

.. _databaseformat:

Database format
-------------------

The first few lines of the IBTrACS database are shown below. ::

    IBTrACS WMO: International Best Tracks Archive for Climate Stewardship -- WMO DATA ONLY -- Version: v03r05
    Serial_Num,Season,Num,Basin,Sub_basin,Name,ISO_time,Nature,Latitude,Longitude,Wind(WMO),Pres(WMO),Center,Wind(WMO) Percentile,Pres(WMO) Percentile,Track_type
    N/A,Year,#,BB,BB,N/A,YYYY-MM-DD HH:MM:SS,N/A,deg_north,deg_east,kt,mb,N/A,%,%,N/A
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-11 06:00:00, NR, -8.60,  79.80,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-12 06:00:00, NR, -9.00,  78.90,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-13 06:00:00, NR,-10.40,  73.20,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-14 06:00:00, NR,-12.80,  69.90,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-15 06:00:00, NR,-13.90,  68.90,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-16 06:00:00, NR,-15.30,  67.70,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-17 06:00:00, NR,-16.50,  67.00,  0.0,    0.0,reunion,-100.000,-100.000,main
    1848011S09080,1848,02, SI, MM,XXXX848003,1848-01-18 06:00:00, NR,-18.00,  67.40,  0.0,    0.0,reunion,-100.000,-100.000,main

.. _configuringsource:

Setting the configuration options
---------------------------------

The columns arrangement can be determined from looking at the second
line of the database:
``Serial_Num,Season,Num,Basin,Sub_basin,Name,ISO_time,Nature,Latitude,Longitude,Wind(WMO),Pres(WMO),Center,Wind(WMO)
Percentile,Pres(WMO) Percentile,Track_type``. The important variables
required for TCRM are the ``Serial_Num``, ``ISO_time``, ``Latitude``,
``Longitude`` and ``Pres(WMO)`` columns. All other variables can be
ignored. Based on this, the ``Columns`` option can be set as::

    Columns = tcserialno,season,num,skip,skip,skip,date,skip,lat,lon,skip,pressure

Notice that there are only 12 columns specified, but the database
contains 15 columns. The last three columns are not used, so since
they are not specified, they are automatically ignored.

The data are comma-delimited, so the ``FieldDelimiter`` option is set
to a comma -- it does not need to be wrapped in quotes. ::
    
    FieldDelimiter = ,


The first three lines are metadata and information that describe the
data (version, column names, units). These lines are ignored when the
``NumberOfHeadingLines`` option is set to 3. ::

    NumberOfHeadingLines = 3

The format of the date field needs to be specified in the
configuration file. Format specification is controlled by Python's
standard :mod:`datetime` module. ::

    DateFormat = %Y-%m-%d %H:%M:%S


The remaining three options -- ``PressureUnits``, ``LengthUnits`` and
``SpeedUnits`` -- set the units of the pressure values, any length
units and the speed units (of motion of the storms). Notice that the
third line in the database indicates the wind speed is in units of
``kt`` (knots) - this refers to the observed maximum wind speed. The
``SpeedUnits`` option controls the units for the forward motion speed
of storms - a field that is not included in the IBTrACS database. ::

    PressureUnits = hPa
    LengthUnits = km
    SpeedUnits = kmh


The full section for using IBTrACS is shown below. ::

    [IBTRACS]
    ; Input data file settings
    URL = ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r05/wmo/csv/Allstorms.ibtracs_wmo.v03r05.csv.gz
    Path = input
    Filename = Allstorms.ibtracs_wmo.v03r05.csv
    NumberOfHeadingLines = 3
    Columns = tcserialno,season,num,skip,skip,skip,date,skip,lat,lon,skip,pressure
    FieldDelimiter = ,
    DateFormat = %Y-%m-%d %H:%M:%S
    PressureUnits = hPa
    LengthUnits = km
    SpeedUnits = kph
