import time
import logging as log
from functools import wraps, reduce

import numpy as np

def fromrecords(records, names):
    """ Convert records to array, even if no data """
    # May become redundant after https://github.com/numpy/numpy/issues/1862
    if records:
        rval = np.rec.fromrecords(records, names=names)
    else:
        rval = np.array([], [(name, 'O') for name in names.split(',')])

    return rval

def timer(func):
    """
    A simple timing decorator for the entire process.

    """

    @wraps(func)
    def wrap(*args, **kwargs):
        t1 = time.time()
        res = func(*args, **kwargs)
        tottime = time.time() - t1
        msg = "%02d:%02d:%02d " % \
            reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
                   [(tottime,), 60, 60])
        log.debug("Time for {0}: {1}".format(func.__name__, msg))
        return res

    return wrap

@timer
def locationRecordsExceeding(hazard_db, locId, windSpeed):
    """
    Select all records where the wind speed at the given location is
    greater than some threshold.

    :param hazard_db: :class:`HazardDatabase` instance.
    :param int locId: Location identifier.
    :param float windSpeed: Select all records where the wind speed
                            at the given location is greater than
                            this value.

    :returns: :class:`numpy.recarray` containing the name, longitude
              & latitude of the location, the wind speed of the
              record, the event Id and the event file that holds the
              event that generated the wind speed.

    Example::

        >>> db = HazardDatabase(configFile)
        >>> locId = 00001
        >>> records = locationRecordsExceeding(db, locId, 47.)

    """

    query = ("SELECT l.locId, l.locName, w.wspd, w.eventId "
             "FROM tblLocations l "
             "INNER JOIN tblWindSpeed w ON l.locId = w.locId "
             "WHERE w.wspd > ? and l.locId = ? "
             "ORDER BY w.wspd ASC")

    cur = hazard_db.execute(query, (windSpeed, locId,))
    results = cur.fetchall()
    results = fromrecords(results, names=('locId,locName,wspd,eventId'))

    return results

@timer
def locationRecords(hazard_db, locId):
    """
    Select all wind speed records for a given location.

    :param hazard_db: :class:`HazardDatabase` instance.
    :param int locId: Location identifier.

    :returns: :class:`numpy.recarray` containing the location id, location
              name, wind speed and event id.

    """

    query = ("SELECT w.locId, l.locName, w.wspd, w.umax, w.vmax, w.eventId "
             "FROM tblWindSpeed w "
             "INNER JOIN tblLocations l "
             "ON w.locId = l.locId "
             "WHERE l.locId = ? ORDER BY w.wspd ASC")
    cur = hazard_db.execute(query, (locId,))
    results = cur.fetchall()
    results = fromrecords(results,
                          names=('locId,locName,wspd,umax,vmax,eventId'))

    return results

@timer
def locationPassage(hazard_db, locId, distance=50):
    """
    Select all records from tblTracks that pass within a defined
    distance of the given location

    :param hazard_db: :class:`HazardDatabase` instance.
    :param int locId: Location identifier.
    :param distance: Distance threshold (in kilometres).

    :returns: :class:`numpy.recarray` containing the location id, location
              name, event id, closest distance of approach, wind speed and
              event file for all events that pass within the defined
              distance of the selected location.

    Example::

        >>> db = HazardDatabase(configFile)
        >>> locId = 000001
        >>> records = locationPassage(db, locId, 50)

    """

    query = ("SELECT l.locId, l.locName, t.eventId, t.distClosest, "
             "w.wspd, e.eventFile FROM tblLocations l "
             "INNER JOIN tblTracks t "
             "ON l.locId = t.locId "
             "JOIN tblWindSpeed w on w.eventId = t.eventId "
             "JOIN tblEvents e on e.eventId = t.eventId "
             "WHERE t.distClosest < ? and l.locId = ?")
    cur = hazard_db.execute(query, (distance, locId))
    results = cur.fetchall()
    results = fromrecords(results,
                          names=('locId,locName,eventId,'
                                 'distClosest,wspd,eventFile'))
    return results

@timer
def locationPassageWindSpeed(hazard_db, locId, speed, distance):
    """
    Select records from _tblWindSpeed_, _tblTracks_ and _tblEvents_ that
    generate a defined wind speed and pass within a given distance
    of the location.

    :param hazard_db: :class:`HazardDatabase` instance.
    :param int locId: Location identifier.
    :param float speed: Minimum wind speed (m/s).
    :param float distance: Distance threshold (kilometres).

    """

    query = ("SELECT l.locName, w.wspd, w.umax, w.vmax, w.eventId, "
             "t.distClosest, e.eventMaxWind, e.eventMinPressure "
             "FROM tblLocations l "
             "JOIN tblWindSpeed w on l.locId = w.locId "
             "JOIN tblEvents e ON e.eventId = w.eventId "
             "JOIN tblTracks t ON w.locId = t.locId AND w.eventId = t.eventId "
             "WHERE l.locId = ? and w.wspd > ? AND t.distClosest <= ? "
             "ORDER BY w.wspd ASC")

    cur = hazard_db.execute(query, (locId, speed, distance))
    results = cur.fetchall()
    results = fromrecords(results,
                          names=('locName,wspd,umax,vmax,eventId,'
                                 'distClosest,maxwind,pmin'))

    return results

@timer
def locationReturnPeriodEvents(hazard_db, locId, return_period):
    """
    Select all records from tblEvents where the wind speed is
    greater than the return period wind speed for the given return period.

    :param hazard_db: :class:`HazardDatabase` instance.
    :param int locId: Location identifier.
    :param int return_period: Nominated return period.

    :returns: :class:`numpy.recarray` of location id and wind speeds of
              all events that are greater than the return level of the
              nominated return period.

    The following example would return the wind speeds of all events that
    exceed the 500-year return period wind speed for the selected location.

    Example::

        >>> db = HazardDatabase(configFile)
        >>> locId = 000001
        >>> records = locationReturnPeriodEvents(db, locId, 500)

    """

    query = ("SELECT l.locId, h.wspd FROM tblLocations l "
             "INNER JOIN tblHazard h ON l.locId = h.locId "
             "WHERE h.returnPeriod = ? and l.locId = ?")
    cur = hazard_db.execute(query, (return_period, locId))
    row = cur.fetchall()
    return_level = row[0][1]
    results = locationRecordsExceeding(hazard_db, locId, return_level)

    return results

@timer
def locationAllReturnLevels(hazard_db, locId):
    """
    Select all return level wind speeds (including upper and lower
    confidence intervals) for a selected location.

    :param hazard_db: :class:`HazardDatabase` instance.
    :param int locId: Location identifier.

    :returns: :class:`numpy.recarray` containing the location id, location
              name, return period, return period windspeed and lower/upper
              estimates of the return period wind speed.

    """

    query = ("SELECT l.locId, l.locName, h.returnPeriod, h.wspd, "
             "h.wspdLower, h.wspdUpper "
             "FROM tblLocations l INNER JOIN tblHazard h "
             "ON l.locId = h.locId "
             "WHERE l.locId = ? "
             "ORDER BY h.returnPeriod")

    cur = hazard_db.execute(query, (locId,))
    results = cur.fetchall()
    results = fromrecords(results,
                          names=('locId,locName,returnPeriod,'
                                 'wspd,wspdLower,wspdUpper'))

    return results

@timer
def selectEvents(hazard_db):
    """
    Select all events from _tblEvents_.

    :param hazard_db: :class:`HazardDatabase` instance.

    :returns: :class:`numpy.recarray` containing the full listing of each
              event in the table.

    """

    query = "SELECT * FROM tblEvents ORDER BY eventMaxWind ASC"
    cur = hazard_db.execute(query)
    results = cur.fetchall()
    names = ("eventNum,eventId,eventFile,eventTrackFile,eventMaxWind,"
             "eventMinPressure,dtTrackFile,dtWindfieldFile,tcrmVer,"
             "Comments,dtCreated")
    results = fromrecords(results, names=names)
    return results
