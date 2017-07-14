===============
Hazard database
===============

To make it easier to acess hazard data, and to enable users to link an
event to a return period wind speed, TCRM stores information on the
individual events, the wind speeds they generate at a set of
locations, and the return levels calculated in the :mod:`hazard`
module.

The database is built on the :mod:`sqlite3` module that is part of the
standard Python library. Extension to a spatial database system has
been identified as an enhancement for future development.

Tables
======

Following is a description of the tables contained in the
database. Tables are created when the database is first
generated. Tables are populated with information from the
simulation. Presently, the tables are not updateable - that is, the
records cannot be modified once inserted into teh database.

Locations table
---------------

The locations table (*tblLocations*) stores the details of the
locations within the model domain. This table has been derived from a
(global) listing of WMO-identified weather stations.

.. tabularcolumns:: |l|c|p{5cm}|
+-------------+----------------------------------+------------+
| Field name  | Field description                | Field type |
+=============+==================================+============+
| locId       | Unique identifier for the        | Integer    |
|             | location.                        |            |
+-------------+----------------------------------+------------+
| locName     | Name of the location.            | Text       |
+-------------+----------------------------------+------------+
| locLon      | Longitude of the location, given | Real       |
|             | in geographic coordinates.       |            |
+-------------+----------------------------------+------------+
| locLat      | Latitude of the location, given  | Real       |
|             | in geographic coordinates.       |            |
+-------------+----------------------------------+------------+
| locElev     | Elevation of the location        | Real       |
|             | (in metres)                      |            |
+-------------+----------------------------------+------------+
| locCountry  | Country of the location          | Text       |
+-------------+----------------------------------+------------+



Events table
------------

The events table (*tblEvents*) stores details of the individual TC
events, such as maximum wind speed and minimum central pressure, as
well as details on the event identifier and files associated with the
event.


.. tabularcolumns:: |l|c|p{5cm}|
+------------------+---------------------------------------------------------------------+------------+
| Field name       | Field description                                                   | Field type |
+==================+=====================================================================+============+
| eventNumber      | Unique identifying number for the event                             |   Integer  |
+------------------+---------------------------------------------------------------------+------------+
| eventId          | Unique identifier for the event                                     |    Text    |
+------------------+---------------------------------------------------------------------+------------+
| eventFile        | File name where the event data (wind and pressure fields) is stored |    Text    |
+------------------+---------------------------------------------------------------------+------------+
| eventTrackFile   | File name where the track data for the event is stored              |    Text    |
+------------------+---------------------------------------------------------------------+------------+
| eventMaxWind     | Maximum sustained wind speed of the event (m/s)                     |    Real    |
+------------------+---------------------------------------------------------------------+------------+
| eventMinPressure | Minimum central pressure of the event (hPa)                         |    Real    |
+------------------+---------------------------------------------------------------------+------------+
| dtTrackFile      | Creation date/time of the track file                                |  Timestamp |
+------------------+---------------------------------------------------------------------+------------+
| dtWindfieldFile  | Creation date/time of the wind field file                           |  Timestamp |
+------------------+---------------------------------------------------------------------+------------+
| tcrmVersion      | TCRM version (usually a Git hash representing the commit)           |    Text    |
+------------------+---------------------------------------------------------------------+------------+
| Comments         | General comments field                                              |    Text    |
+------------------+---------------------------------------------------------------------+------------+
| dtCreated        | Creation date/time for the record                                   |  Timestamp |
+------------------+---------------------------------------------------------------------+------------+

Wind speed table
----------------

The wind speed table (*tblWindSpeed*) contains records for each
location and event, recording the maximum wind speed at the location,
the component wind speeds (east/west and north/south) and the minimum
pressure at the location.


.. tabularcolumns:: |l|c|p{5cm}|
+------------+--------------------------------------------------------------------+------------+
| Field name | Field description                                                  | Field type |
+============+====================================================================+============+
| locId      | Unique identifier for location                                     | Integer    |
+------------+--------------------------------------------------------------------+------------+
| eventId    | Unique identifier for event                                        | Text       |
+------------+--------------------------------------------------------------------+------------+
| wspd       | Maximum wind speed from the event at the location (m/s)            | Real       |
+------------+--------------------------------------------------------------------+------------+
| umax       | Eastward component of the maximum wind speed from the event (m/s)  | Real       |
+------------+--------------------------------------------------------------------+------------+
| vmax       | Northward component of the maximum wind speed from the event (m/s) | Real       |
+------------+--------------------------------------------------------------------+------------+
| pmin       | Minimum pressure from the event at the location (hPa)              | Real       |
+------------+--------------------------------------------------------------------+------------+
| Comments   | General comments field                                             | Text       |
+------------+--------------------------------------------------------------------+------------+
| dtCreated  | Creation date/time for the record                                  | Timestamp  |
+------------+--------------------------------------------------------------------+------------+


Hazard table
------------

The hazard table (*tblHazard*) records the return period wind speeds
at each location, as well as upper/lower estimated confidence
intervals (if calculated) and the location, scale and shape parameters
of the fitted GEV distribution.

.. tabularcolumns:: |l|c|p{5cm}|
+--------------+--------------------------------------------------------------------------------------+------------+
| Field name   | Field description                                                                    | Field type |
+==============+======================================================================================+============+
| locId        | Unique identifier for location                                                       | Integer    |
+--------------+--------------------------------------------------------------------------------------+------------+
| returnPeriod | Return period (years)                                                                | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| wspd         | Return period wind speed for the location and given return period (m/s)              | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| wspdUpper    | Estimated upper confidence bound of return period wind speed (95th percentile) (m/s) | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| wspdLower    | Estimated lower confidence bound of return period wind speed (5th percentile) (m/s)  | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| loc          | Location parameter for the fitted Generalised Extreme Value (GEV) distribution       | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| scale        | Scale parameter for the fitted GEV distribution                                      | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| shape        | Shape parameter for the fitted GEV distribution                                      | Real       |
+--------------+--------------------------------------------------------------------------------------+------------+
| tcrmVersion  | TCRM version (usually a Git hash representing the commit)                            | Text       |
+--------------+--------------------------------------------------------------------------------------+------------+
| dtHazardFile | Creation date/time for the hazard file                                               | Timestamp  |
+--------------+--------------------------------------------------------------------------------------+------------+
| Comments     | General comments field                                                               | Text       |
+--------------+--------------------------------------------------------------------------------------+------------+
| dtCreated    | Creation date/time for the record                                                    | Timestamp  |
+--------------+--------------------------------------------------------------------------------------+------------+


Tracks table
------------

The tracks table (*tblTracks*) records information on the individual
tracks and their proximity to each individual location in the domain.

.. tabularcolumns:: |l|c|p{5cm}|
+-------------+----------------------------------------------------------------------+------------+
| Field name  | Field description                                                    | Field type |
+=============+======================================================================+============+
| locId       | Unique identifier for location                                       | Integer    |
+-------------+----------------------------------------------------------------------+------------+
| eventId     | Unique identifier for the event                                      | Text       |
+-------------+----------------------------------------------------------------------+------------+
| distClosest | Distance of closest approach of the track to the location (km)       | Real       |
+-------------+----------------------------------------------------------------------+------------+
| prsClosest  | Central pressure of the event at the point of closest approach (hPa) | Real       |
+-------------+----------------------------------------------------------------------+------------+
| dtClosest   | Date/time of the point of closest approach                           | Timestamp  |
+-------------+----------------------------------------------------------------------+------------+
| Comments    | General comments field                                               | Text       |
+-------------+----------------------------------------------------------------------+------------+
| dtCreated   | Creation date/time for the record                                    | Timestamp  |
+-------------+----------------------------------------------------------------------+------------+


Queries
=======

A small number of queries are pre-built into the module. These are
provided to indicate approaches users can use to develop their own
queries that satisfy their needs.

Queries use SQL syntax, and are passed to the
:class:`HazardDatabase.execute` method to perform the SQL command. Queries
are assembled using the :class:`sqlite3` DB-API parameter substitution
(see the `sqlite3 <https://docs.python.org/2/library/sqlite3.html>`_
documentation for more details).

Results from queries are returned as a :class:`numpy.recarray`, with
the dtype defined by the corresponding type of the database field (see
above for field types).

