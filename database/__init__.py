"""
:mod:`database` -- build, update and query a database of hazard and events
==========================================================================

.. module:: database
    :synopsis: Build, update, query a database for hazard and event
               information.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

The :class:`HazardDatabase` class provides methods to create and
populate a :mod:`sqlite3` database that holds location-specific
information on the synthetic events and the hazard (return period wind
speeds). This database allows users to identify events that generate
wind speeds corresponding to some threshold, such as a return period
wind speed, at each location in the domain, which could then be
selected for more detailed modelling.

:class:`HazardDatabase` is initially intended to be created once, then
queried from subsequent scripts or interactive sessions. By
inheriting from the :class:`Singleton` class, this will only ever return
the first instance of :class:`HazardDatabase` created in a session.


TODO::

- Upgrade to spatial database to better handle geometries & projected
  location data.
- Check for existing database and create new/replace existing db if
  configuration settings have changed. Requires config settings to
  be stored in the db (?).
- Separate the query functions to separate files.

"""

import os
import logging
import sqlite3
from sqlite3 import PARSE_DECLTYPES, PARSE_COLNAMES

from datetime import datetime
import unicodedata
import re
from os.path import join as pjoin

from shapely.geometry import Point
from netCDF4 import Dataset
import numpy as np

from Utilities.config import ConfigParser
from Utilities.files import flGetStat
from Utilities.maputils import find_index
from Utilities.track import loadTracksFromFiles, loadTracks
from Utilities.singleton import Singleton
from Utilities.parallel import attemptParallel, disableOnWorkers
from Utilities.process import pAlreadyProcessed, pWriteProcessedFile, \
    pGetProcessedFiles
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

# pylint: disable=R0914,R0902

# Table definition statements
# Stations - we assume a geographic coordinate system:
TBLLOCATIONDEF = ("CREATE TABLE IF NOT EXISTS tblLocations "
                  "(locId integer PRIMARY KEY, locCode text, "
                  "locName text, locType text, locLon real, "
                  "locLat real, locElev real, locCountry text, "
                  "locSource text, Comments text, "
                  "dtCreated timestamp)")

# Events:
TBLEVENTSDEF = ("CREATE TABLE IF NOT EXISTS tblEvents "
                "(eventNumber integer PRIMARY KEY, eventId text, "
                "eventFile text, eventTrackFile text, "
                "eventMaxWind real, eventMinPressure real, "
                "dtTrackFile timestamp, dtWindfieldFile timestamp, "
                "tcrmVersion text, Comments text, dtCreated timestamp)")

#Station wind speed from events:
TBLWINDSPEEDDEF = ("CREATE TABLE IF NOT EXISTS tblWindSpeed "
                   "(locId integer, eventId text, wspd real, umax real, "
                   "vmax real, pmin real, Comments text, "
                   "dtCreated timestamp)")

# Station hazard levels:
TBLHAZARDDEF = ("CREATE TABLE IF NOT EXISTS tblHazard "
                "(locId integer, returnPeriod real, wspd real, "
                " wspdUpper real, wspdLower real, loc real, "
                "scale real, shape real, tcrmVersion text, "
                "dtHazardFile timestamp, Comments text, "
                "dtCreated timestamp)")

# Proximity of tracks to stations:
TBLTRACKSDEF = ("CREATE TABLE IF NOT EXISTS tblTracks "
                "(locId integer, eventId text, distClosest real, "
                "prsClosest real, dtClosest timestamp, Comments text, "
                "dtCreated timestamp)")

# Insert statements:
# Insert locations:
INSLOCATIONS = ("INSERT OR REPLACE INTO tblLocations "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?)")

# Insert event record:
INSEVENTS = "INSERT INTO tblEvents VALUES (?,?,?,?,?,?,?,?,?,?,?)"

# Insert wind speed record:
INSWINDSPEED = ("INSERT INTO tblWindSpeed (locId, eventId, wspd, "
                "umax, vmax, pmin, Comments, dtCreated) "
                "VALUES (?,?,?,?,?,?,?,?)")

# Insert hazard record:
INSHAZARD = "INSERT INTO tblHazard VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"

# Insert track record:
INSTRACK = "INSERT INTO tblTracks VALUES (?,?,?,?,?,?,?)"

# Select statements;
# Select locations within domain:
SELECTLOCATIONS = ("SELECT * FROM tblLocations WHERE "
                   "locLon >= ? and locLon <= ? and "
                   "locLat >= ? and locLat <= ?")

# Select locId, locLon & locLat from the subset of locations:
SELECTLOCLONLAT = "SELECT locId, locLon, locLat FROM tblLocations "

def windfieldAttributes(ncobj):
    """
    Extract the required attributes from a netCDF file.

    :param str ncobj: :class:`netCDF4.Dataset` instance.

    :returns: A tuple containing the track filename, file modification date,
              TCRM version, minimum pressure and maximum wind, stored as
              global attributes in the netCDF file.

    """

    trackfile = getattr(ncobj, 'track_file')
    trackfiledate = getattr(ncobj, 'track_file_date')
    trackfiledate = datetime.strptime(trackfiledate, '%Y-%m-%d %H:%M:%S')
    tcrm_vers = getattr(ncobj, 'tcrm_version')
    trackfile = unicodedata.normalize("NFKD",
                                      trackfile).encode('utf-8', 'ignore')
    trackfile = os.path.basename(trackfile)
    tcrm_vers = unicodedata.normalize("NFKD",
                                      tcrm_vers).encode('utf-8', 'ignore')
    slpobj = ncobj.variables['slp']
    minslp = getattr(slpobj, 'actual_range')[0]

    vmaxobj = ncobj.variables['vmax']
    maxwind = getattr(vmaxobj, 'actual_range')[1]

    return (trackfile, trackfiledate, tcrm_vers, minslp, maxwind)

def HazardDatabase(configFile): # pylint: disable=C0103
    """
    A wrapper function to instantiate :class:`_HazardDatabase`. We actually
    call :meth:`getInstance` to get a singleton instance of the
    :class:`_HazardDatabase`. If one exists already, then that instance is
    returned.

    :param str configFile: Path to configuration file

    """
    return _HazardDatabase.getInstance(configFile)

class _HazardDatabase(sqlite3.Connection, Singleton):
    """
    Create and update a database of locations, events, hazard, wind
    speed and tracks. Because it subclasses the :class:`Singleton` object,
    it will create a new instance if there is not already an active instance,
    or it will return an existing instance.

    :param str configFile: Path to the simulation configuration file.

    """
    ignoreSubsequent = True
    def __init__(self, configFile):

        config = ConfigParser()
        config.read(configFile)

        self.outputPath = config.get('Output', 'Path')
        self.windfieldPath = pjoin(self.outputPath, 'windfield')
        self.trackPath = pjoin(self.outputPath, 'tracks')
        self.hazardPath = pjoin(self.outputPath, 'hazard')
        self.domain = config.geteval('Region', 'gridLimit')
        self.hazardDB = pjoin(self.outputPath, 'hazard.db')
        self.locationDB = pjoin(self.outputPath, 'locations.db')
        self.datfile = config.get('Process', 'DatFile')
        self.excludePastProcessed = config.getboolean('Process',
                                                      'ExcludePastProcessed')

        pGetProcessedFiles(self.datfile)

        sqlite3.Connection.__init__(self, self.hazardDB,
                                    detect_types=PARSE_DECLTYPES|PARSE_COLNAMES)

        self.exists = True

        import atexit
        atexit.register(self.close)

    @disableOnWorkers
    def createDatabase(self):
        """
        Create the database and the tables in the database.

        """
        log.info("Building the hazard database...")
        self.createTable('tblLocations', TBLLOCATIONDEF)
        self.createTable('tblEvents', TBLEVENTSDEF)
        self.createTable('tblWindSpeed', TBLWINDSPEEDDEF)
        self.createTable('tblHazard', TBLHAZARDDEF)
        self.createTable('tblTracks', TBLTRACKSDEF)
        self.exists = True
        self.commit()
        return

    @disableOnWorkers
    def createTable(self, tblName, tblDef):
        """
        Create a table.

        :param tblName: Table name.
        :param tblDef: Table definition.

        :raises: `sqlite3.Error` if unable to create the database.
        """
        log.info("Creating table {0}".format(tblName))
        log.debug("Executing statement: {0}".format(tblDef))
        try:
            self.execute(tblDef)
            self.commit()
        except sqlite3.Error as err:
            log.exception("Cannot create table {0}: {1}".\
                          format(tblName, err.args[0]))
            raise

    @disableOnWorkers
    def setLocations(self):
        """
        Populate _tblLocations_ in the hazard database with all
        locations from the default locations database that lie
        within the simulation domain. If the table exists and is
        populated, the records will be updated and any new records
        inserted.

        """

        conn = sqlite3.connect(self.locationDB,
                               detect_types=PARSE_DECLTYPES|PARSE_COLNAMES)
        cur = conn.execute(SELECTLOCATIONS, (self.domain['xMin'],
                                             self.domain['xMax'],
                                             self.domain['yMin'],
                                             self.domain['yMax']))

        locations = cur.fetchall()
        conn.close()
        if len(locations) >= 1:
            self.executemany(INSLOCATIONS, locations)
            self.commit()
        else:
            log.info("No locations returned")

    def getLocations(self):
        """
        Retrieve all locations stored in the hazard database.

        :returns: List of tuples containing location id, longitude and latitude.

        :raises: `sqlite3.Error` if unable to retrieve the locations.
        """
        try:
            cur = self.execute(("SELECT locId, locName, locLon, locLat "
                                "FROM tblLocations"))
        except sqlite3.Error as err:
            log.exception("Cannot retrieve locations from tblLocations: {0}".\
                          format(err.args[0]))
            raise
        else:
            locations = cur.fetchall()

        locations = np.rec.fromrecords(locations,
                                       names=("locId,locName,locLon,locLat"))
        return locations

    def generateEventTable(self):
        """
        Populate _tblEvents_ with the details of the synthetic events
        generated in the simulation. This table only holds the
        metadata of the events. At this time, since TCRM generates
        annual event sets, this table only stores details of the
        annual event set. Future versions can hold additional metadata
        about each individual synthetic TC event.

        :raises: `sqlite3.Error` if unable to insert events into
                 _tblEvents_.

        """
        log.info("Inserting records into tblEvents")

        fileList = os.listdir(self.windfieldPath)
        fileList = [f for f in fileList if
                    os.path.isfile(pjoin(self.windfieldPath, f))]

        pattern = re.compile(r'\d+')
        params = []
        for n, f in enumerate(sorted(fileList)):
            log.debug("Processing {0}".format(f))
            sim, num = pattern.findall(f)
            eventId = "%03d-%05d" % (int(sim), int(num))
            fname = pjoin(self.windfieldPath, f)
            si = os.stat(fname)
            dtWindfieldFile = datetime.fromtimestamp(int(si.st_mtime))
            trackfile, dtTrackFile, tcrm_version, minslp, maxwind = \
                windfieldAttributes(fname)
            params.append(("%06d"%n, eventId, os.path.basename(fname),
                           trackfile, float(maxwind), float(minslp),
                           dtTrackFile, dtWindfieldFile, tcrm_version,
                           "", datetime.now()))


        try:
            self.executemany(INSEVENTS, params)
        except sqlite3.Error as err:
            log.exception("Cannot insert records into tblEvents: {0}".\
                          format(err.args[0]))
            raise
        else:
            self.commit()

    def insertEvents(self, eventparams):
        """
        Insert records into _tblEvents_, using the collection of event
        parameters.

        :param eventparams: A collection of tuples, each of which represents
                            an individual event record.
        """
        log.debug("Inserting records into tblEvents")
        try:
            self.execute(INSEVENTS, eventparams)
        except sqlite3.Error as err:
            log.exception("Cannot insert records into tblEvents: {0}".\
                          format(err.args[0]))
        else:
            self.commit()

    def insertWindSpeeds(self, wsparams):
        """
        Insert records into _tblWindSpeed_, using the collection of
        location-specific wind speeds extracted from the event catalogue.

        :param wsparams: A collection of tuples, each of which represents
                         an individual location record.
        """
        nrecords = len(wsparams)
        log.debug("Inserting {0} records into tblWindSpeed".format(nrecords))
        try:
            self.executemany(INSWINDSPEED, wsparams)
        except sqlite3.Error as err:
            log.exception("Cannot insert records into tblWindSpeed: {0}".\
                          format(err.args[0]))
        except sqlite3.ProgrammingError as err:
            log.exception("Programming error: {0}".format(err.args[0]))
        else:
            self.commit()


    def processEvents(self):
        """
        Process the events (wind fields) for each location within the
        model domain and populate _tblWindSpeed_. This will store the
        modelled wind speed (or the missing value) at each grid point,
        from each synthetic event.

        :raises: `sqlite3.Error` if unable to insert records into
                 _tblWindSpeed_.

        """
        fileList = os.listdir(self.windfieldPath)
        fileList = [f for f in fileList if
                    os.path.isfile(pjoin(self.windfieldPath, f))]

        work_tag = 0
        result_tag = 1
        if (pp.rank() == 0) and (pp.size() > 1):
            locations = self.getLocations()
            w = 0
            p = pp.size() - 1
            for d in range(1, pp.size()):
                if w < len(fileList):
                    pp.send((fileList[w], locations, w),
                            destination=d, tag=work_tag)
                    log.debug("Processing {0} ({1} of {2})".\
                              format(fileList[w], w, len(fileList)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    p = w

            terminated = 0

            while terminated < p:
                result, status = pp.receive(pp.any_source, tag=result_tag,
                                            return_status=True)
                log.debug("Processing results from node {0}".\
                          format(status.source))
                eventparams, wsparams = result
                self.insertEvents(eventparams)
                self.insertWindSpeeds(wsparams)
                log.debug("Done inserting records from node {0}".\
                          format(status.source))

                d = status.source
                if w < len(fileList):
                    pp.send((fileList[w], locations, w),
                            destination=d, tag=work_tag)
                    log.debug("Processing file %d of %d" % (w, len(fileList)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    terminated += 1

        elif (pp.size() > 1) and (pp.rank() != 0):
            while True:
                W = pp.receive(source=0, tag=work_tag)
                if W is None:
                    break

                log.info("Processing {0} on node {1}".\
                         format(W[0], pp.rank()))
                results = self.processEvent(*W)
                log.debug("Results received on node {0}".format(pp.rank()))
                pp.send(results, destination=0, tag=result_tag)

        elif pp.size() == 1 and pp.rank() == 0:
            # Assume no Pypar:
            locations = self.getLocations()
            for eventNum, filename in enumerate(fileList):
                log.debug("Processing {0} ({1} of {2})".format(filename, 
                                                               eventNum, 
                                                               len(fileList)))
                result = self.processEvent(filename, locations, eventNum)
                eventparams, wsparams = result
                self.insertEvents(eventparams)
                self.insertWindSpeeds(wsparams)


    def loadWindfieldFile(self, ncobj):
        """
        Load an individual dataset.

        :param str filename: filename to load.

        :returns: tuple containing longitude, latitude, wind speed,
                  eastward and northward components and pressure grids.
        """
        lon = ncobj.variables['lon'][:]
        lat = ncobj.variables['lat'][:]
        vmax = ncobj.variables['vmax'][:]
        ua = ncobj.variables['ua'][:]
        va = ncobj.variables['va'][:]
        pmin = ncobj.variables['slp'][:]

        return (lon, lat, vmax, ua, va, pmin)

    def processEvent(self, filename, locations, eventNum):
        """
        Process an individual event file
        :param str filename: Name of a file to process.
        :param list locations: List of locations to sample data for.
        :param int eventNum: Ordered event number.
        """
        log.debug("Processing {0}".format(pjoin(self.windfieldPath, filename)))
        pattern = re.compile(r'\d+')
        sim, num = pattern.findall(filename)
        eventId = "%03d-%05d" % (int(sim), int(num))
        log.debug("Event ID: {0}".format(eventId))
        try:
            ncobj = Dataset(pjoin(self.windfieldPath, filename))
        except:
            log.warn("Cannot open {0}".\
                     format(pjoin(self.windfieldPath, filename)))

        # First perform the event update for tblEvents:
        fname = pjoin(self.windfieldPath, filename)
        log.debug("Filename: {0}".format(fname))
        si = os.stat(fname)
        dtWindfieldFile = datetime.fromtimestamp(int(si.st_mtime))
        trackfile, dtTrackFile, tcrm_version, minslp, maxwind = \
                windfieldAttributes(ncobj)

        eventparams = ("%06d"%eventNum, eventId, os.path.basename(fname),
                       trackfile, float(maxwind), float(minslp),
                       dtTrackFile, dtWindfieldFile, tcrm_version,
                       "", datetime.now())

        # Perform update for tblWindSpeed:
        lon, lat, vmax, ua, va, pmin = self.loadWindfieldFile(ncobj)
        ncobj.close()
        wsparams = list()

        for loc in locations:
            locId, locName, locLon, locLat = loc
            i = find_index(lon, locLon)
            j = find_index(lat, locLat)
            locVm = vmax[j, i]
            locUa = ua[j, i]
            locVa = va[j, i]
            locPr = pmin[j, i]
            locParams = (locId, eventId, float(locVm), float(locUa),
                         float(locVa), float(locPr), " ", datetime.now())
            wsparams.append(locParams)

        log.debug("Finished extracting data from {0}".format(filename))
        return (eventparams, wsparams,)

    @disableOnWorkers
    def processHazard(self):
        """
        Update _tblHazard_ with the return period wind speed data.

        """
        log.info("Inserting records into tblHazard")
        locations = self.getLocations()
        hazardFile = pjoin(self.hazardPath, 'hazard.nc')
        ncobj = Dataset(hazardFile)

        try:
            tcrm_version = getattr(ncobj, 'tcrm_version')
        except AttributeError:
            log.info("Missing tcrm_version attribute from {0}".\
                     format(hazardFile))
            tcrm_version = ''

        si = os.stat(hazardFile)
        dtHazardFile = datetime.fromtimestamp(int(si.st_mtime))
        lon = ncobj.variables['lon'][:]
        lat = ncobj.variables['lat'][:]
        years = ncobj.variables['years'][:]

        wspd = ncobj.variables['wspd'][:]
        wspdUpper = ncobj.variables['wspdupper'][:]
        wspdLower = ncobj.variables['wspdlower'][:]
        locationParam = ncobj.variables['loc'][:]
        scaleParam = ncobj.variables['scale'][:]
        shpParam = ncobj.variables['shp'][:]
        ncobj.close()

        params = []
        for k, year in enumerate(years):
            for loc in locations:
                locId, locName, locLon, locLat = loc
                log.debug("Extracting data for location: {0}".format(locName))
                i = find_index(lon, locLon)
                j = find_index(lat, locLat)
                locWspd = wspd[k, j, i]
                locUpper = wspdUpper[k, j, i]
                locLower = wspdLower[k, j, i]

                locLoc = locationParam[j, i]
                locScale = scaleParam[j, i]
                locShp = shpParam[j, i]

                locParams = (locId, int(year), float(locWspd), float(locUpper),
                             float(locLower), float(locLoc), float(locScale),
                             float(locShp), tcrm_version, dtHazardFile, "",
                             datetime.now())

                params.append(locParams)

        try:
            self.executemany(INSHAZARD, params)
        except sqlite3.Error as err:
            log.exception("Cannot insert records into tblHazard: {0}".\
                          format(err.args[0]))
            raise
        else:
            self.commit()

    def processTracks(self):
        """
        Populate tblTracks with the details of tracks and their proximity to
        the locations in the domain.

        """
        log.info("Inserting records into tblTracks")
        locations = self.getLocations()

        log.debug("Processing tracks in {0}".format(self.trackPath))
        files = os.listdir(self.trackPath)
        trackfiles = [pjoin(self.trackPath, f) for f in files if f.startswith('tracks')]
        #tracks = [t for t in loadTracksFromFiles(sorted(trackfiles))]
        log.debug("There are {0} track files".format(len(trackfiles)))


        work_tag = 0
        result_tag = 1

        if (pp.rank() == 0) and (pp.size() > 1):
            w = 0
            p = pp.size() - 1
            for d in range(1, pp.size()):
                if w < len(trackfiles):
                    pp.send((trackfiles[w], locations),
                            destination=d, tag=work_tag)
                    log.debug("Processing {0}".format(trackfiles[w]))
                    log.debug("Processing track {0:d} of {1:d}".\
                              format(w, len(trackfiles)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    p = w

            terminated = 0

            while terminated < p:
                result, status = pp.receive(pp.any_source, tag=result_tag,
                                            return_status=True)
                if result:
                    log.debug("Inserting results into tblTracks")
                    self.insertTracks(result)
                d = status.source
                if w < len(trackfiles):
                    pp.send((trackfiles[w], locations),
                            destination=d, tag=work_tag)
                    log.debug("Processing track {0:d} of {1:d}".\
                              format(w, len(trackfiles)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    terminated += 1

        elif (pp.size() > 1) and (pp.rank() != 0):
            while True:
                W = pp.receive(source=0, tag=work_tag)
                log.debug("Received track on node {0}".format(pp.rank()))
                if W is None:
                    break

                results = self.processTrack(*W)
                pp.send(results, destination=0, tag=result_tag)

        elif pp.size() == 1 and pp.rank() == 0:
            # No Pypar
            for w, trackfile in enumerate(trackfiles):
                log.debug("Processing trackfile {0:d} of {1:d}".\
                          format(w, len(trackfiles)))
                #if len(track.data) == 0:
                #    continue
                result = self.processTrack(trackfile, locations)
                if result is not None:
                    self.insertTracks(result)

    def processTrack(self, trackfile, locations):
        """
        Process individual track to determine distance to locations, etc.

        Any empty tracks are filtered at an earlier stage. If an empty
        track is passed, then a None result is returned.

        :param track: :class:`Track` instance.
        :param locations: list of locations in the simulation domain.
        """
        tracks = [t for t in loadTracksFromFiles([trackfile])]
        points = [Point(loc[2], loc[3]) for loc in locations]
        records = []
        for track in tracks:
            if len(track.data) == 0:
                log.debug("Got an empty track: returning None")
                continue #return None
            distances = track.minimumDistance(points)
            for (loc, dist) in zip(locations, distances):
                locRecs = (loc[0], "%03d-%05d"%(track.trackId),
                           dist, None, None, "", datetime.now())

                records.append(locRecs)
            log.debug("Track {0}-{1} has {2} records".format(track.trackId[0], 
                                                             track.trackId[1],
                                                             len(records)))
        return records

    def insertTracks(self, trackRecords):
        """
        Insert track parameters into _tblTracks_, using details for each
        location in relation to the track

        :param trackParams: collection of tuples that hold records (i.e.
                            (location details) for each track.
        """
        try:
            self.executemany(INSTRACK, trackRecords)
        except sqlite3.Error as err:
            log.exception(("Cannot insert records into tblTracks: "
                           "{0}").format(err.args[0]))
            raise
        else:
            self.commit()
            log.debug("Inserted {0} records into tblTracks".format(len(trackRecords)))

def run(configFile):
    """
    Run database update

    :param str configFile: path to a configuration file.

    """

    log.info('Running database update')

    config = ConfigParser()
    config.read(configFile)
    outputPath = config.get('Output', 'Path')
    location_db = pjoin(outputPath, 'locations.db')
    if not os.path.exists(location_db):
        location_file = config.get('Input', 'LocationFile')
        buildLocationDatabase(location_db, location_file)

    global pp
    pp = attemptParallel()

    db = HazardDatabase(configFile)

    db.createDatabase()
    db.setLocations()

    pp.barrier()
    db.processEvents()
    pp.barrier()

    db.processHazard()

    pp.barrier()
    db.processTracks()
    pp.barrier()

    db.close()
    pp.barrier()
    log.info("Created and populated database")
    log.info("Finished running database creation")

@disableOnWorkers
def buildLocationDatabase(location_db, location_file, location_type='AWS'):
    """
    Build a database of locations, using a point shape file of the locations.
    The locations *must* be represented in a geographic coordinate system.

    This version is hard coded to work with the `stationlist` file that is
    provided with the RIP4 graphics package, which has in turn been stored
    as a shapefile.

    Users can augment the basic location database with their own data, noting
    the schema for ``tblLocations``.

    :param str location_db: Path to the location database.
    :param str location_file: Path to a shape file containing location data.

    :returns: :class:`numpy.recarray` containing location Id, name, longitude,
              latitude, elevation, country, comments and current datetime.

    TODO: Build a way to ingest user-defined list of fields that correspond
          to the required fields in tblLocations. e.g. using a mappings dict::

              mappings = {
                  'locCode' : 'WMO',
                  'locName' : 'Place'
                  'locCountry' : 'Cou'
                  'locElev' : 'Elevation'
                  'Comments' : 'ICAO'
                  }

              columns = ('locCode', 'locName', 'locCountry',
                         'locElev', 'Comments')

              for col in columns:
                  if mappings.has_key(col):
                      field = shpGetField(location_file, mappings[col])

    """

    from Utilities.shptools import shpReadShapeFile
    log.info("Creating location database")
    locations = []
    vertices, records = shpReadShapeFile(location_file)

    # Perform a check that locations are in geographic coordinates:
    lons = []
    lats = []
    for v in vertices.values():
        lon, lat = v[0]
        lons.append(lon)
        lats.append(lat)

    msg = ("Location shapefile must be in a geograpic coordinate system "
           "(i.e. it must have lat/lon vertices). It looks like this "
           "one has vertices in map projection coordinates. You can convert "
           "the shapefile to geographic coordinates using the shpproj utility "
           "from the shapelib tools "
           "(http://shapelib.maptools.org/shapelib-tools.html)")

    if (max(lons) > 721.) or (min(lons) < -721.) or \
         (max(lats) > 91.) or (min(lats) < -91):
        raise ValueError(msg)

    # Prepare entries:
    for v, r in zip(vertices.values(), records):
        locLon, locLat = v[0]
        locLon = np.mod(locLon, 360.)
        locCode = str(r[0])
        locName = r[2]
        locCountry = r[4]
        locElev = r[7]
        locComment = r[1]
        log.debug("Inserting record for: {0} ({1}): ({2}, {3})".\
                      format(locName, locCode, locLon, locLat))
        locations.append((None, locCode, locName, location_type,
                          locLon, locLat, locElev, locCountry,
                          os.path.basename(location_file),
                          locComment, datetime.now()))

    locdb = sqlite3.connect(location_db,
                            detect_types=PARSE_DECLTYPES|PARSE_COLNAMES)
    locdb.execute(TBLLOCATIONDEF)
    locdb.executemany(INSLOCATIONS, locations)
    locdb.commit()
    locdb.close()

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
    results = np.rec.fromrecords(results,
                                 names=('locId,locName,wspd,eventId'))

    return results

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
    results = np.rec.fromrecords(results,
                                 names=('locId,locName,wspd,'
                                        'umax,vmax,eventId'))

    return results

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
    results = np.rec.fromrecords(results,
                                 names=('locId,locName,eventId,'
                                        'distClosest,wspd,eventFile'))
    return results

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
    results = np.rec.fromrecords(results,
                                 names=('locName,wspd,umax,vmax,eventId,'
                                        'distClosest,maxwind,pmin'))

    return results

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
    results = np.rec.fromrecords(results,
                                 names=('locId,locName,returnPeriod,'
                                        'wspd,wspdLower,wspdUpper'))

    return results

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
    results = np.rec.fromrecords(results, names=names)
    return results
