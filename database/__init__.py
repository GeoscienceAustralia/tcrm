"""
:mod:`database` -- build, update and query a database of hazard and events
==========================================================================

.. module:: database
    :synopsis: Build, update, query a database for hazard and event
               information.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>


"""

import os
import itertools
import logging
import sqlite3
from datetime import datetime
import unicodedata

from os.path import join as pjoin

from netCDF4 import Dataset

from Utilities.config import ConfigParser
from Utilities.files import flModDate
from Utilities.maputils import find_index

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

# Stations - we assume a geographic corrdinate system:
tblLocationsDef = ("CREATE TABLE tblLocations "
                  "(locId text, locName text, locLon real, "
                  "locLat real, locElev real, locCountry text, "
                  "Comments text, dtCreated timestamp)")

# Events:
tblEventsDef = ("CREATE TABLE tblEvents "
                "(eventId text, eventFile text, eventTrackFile text, "
                "eventMaxWind real, eventMinPressure real, "
                "dtTrackFile timestamp, dtWindfieldFile timestamp, "
                "tcrmVersion text, Comments text, dtCreated timestamp)")

#Station wind speed from events:
tblWindSpeedDef = ("CREATE TABLE tblWindSpeed "
                   "(locId text, eventId text, wspd real, umax real, "
                   "vmax real, pmin real, Comments text, "
                   "dtCreated timestamp)")

# Station hazard levels:
tblHazardDef = ("CREATE TABLE tblHazard "
                "(locId text, returnPeriod real, wspd real, "
                " wspdUpper real, wspdLower real, loc real, "
                "scale real, shape real, tcrmVersion text, "
                "dtHazardFile timestamp, Comments text, "
                "dtCreated timestamp)")

# Proximity of tracks to stations:
tblTracksDef = ("CREATE TABLE tblTracks "
                "(locId text, eventId text, distClosest real, "
                "prsClosest real, dtClosest timestamp, Comments text, "
                "dtCreated timestamp)")

# Insert locations:
insLocations = "INSERT INTO tblLocations VALUES (?,?,?,?,?,?,?,?)"

# Insert event record:
insEvents = "INSERT INTO tblEvents VALUES (?,?,?,?,?,?,?,?,?,?)"

# Insert wind speed record:
insWindSpeed = "INSERT INTO tblWindSpeed VALUES (?,?,?,?,?,?,?,?)"

# Insert hazard record:
insHazard = "INSERT INTO tblHazard VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"

# Insert track record:
insTrack = "INSERT INTO tblTracks VALUES (?,?,?,?,?,?)"

# Select locations within domain:
selectLocations = ("SELECT * FROM tblLocations WHERE "
                   "locLon >= ? and locLon <= ? and "
                   "locLat >= ? and locLat <= ?")

# Select locId, locLon & locLat from the subset of locations:
selectLocLonLat = "SELECT locId, locLon, locLat FROM tblLocations "

def singleton(cls):
    instances = {}

    def getinstance(*args, **kwargs):
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance


@singleton
class database(object):
    """
    Create and update a database of locations, events, hazard, wind speed and tracks.

    :param str dbFile: Path to the database file.

    """

    def __init__(self, dbFile):

        self.dbFile = dbFile
        self.exists = False
        if os.path.isfile(self.dbFile):
            self.exists = True

    def close(self):
        self.conn.close()

    def connect(self):
         self.conn = sqlite3.connect(self.dbFile,
                                     detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)

    def create_database(self):
        """
        Create the database.

        """
        if self.exists:
            log.warn("Database file already exists -- passing")
            pass

        else:
            try:
                self.conn = sqlite3.connect(self.dbFile)
            except sqlite3.Error as e:
                log.exception("Cannot create database: %s"% e.args[0])
            else:
                self.exists = True

        if self.exists:
            # Database file exists -- skip creation
            log.debug("Database file already exists -- passing")
            self.conn = sqlite3.connect(self.dbFile)
            return
        else:
            log.info("Building the hazard database...")
            self.conn = sqlite3.connect(self.dbFile)
            self.create_table('tblLocations', tblLocationsDef)
            self.create_table('tblEvents', tblEventsDef)
            self.create_table('tblWindSpeed', tblWindSpeedDef)
            self.create_table('tblHazard', tblHazardDef)
            self.create_table('tblTracks', tblTracksDef)
            self.exists = True
            self.conn.commit()
            return

    def create_table(self, tblName, tblDef):
        """
        Create a table.

        :param tblName: Table name.
        :param tblDef: Table definition.

        """
        log.info("Creating table %s" % tblName)
        log.debug("Executing statement: %s" % tblDef)
        try:
            c = self.conn.cursor()
            c.execute(tblDef)
            self.conn.commit()
        except sqlite3.Error as e:
            log.exception("Cannot create table: %s" % e.args[0])

    def execute(self, statement, parameters=None):
        """
        Update the database using the statement and values.

        :param str statement: SQL statement.
        :param parameters: Tuple of values (for single statement)
                           or list of tuples containing values to
                           update/insert into the database.

        """
        log.debug("Executing statement: %s" % statement)
        cursor = self.conn.cursor()
        if type(parameters) == list and len(parameters > 1):
            # Multiple values:
            cursor.executemany(statement, parameters)
        else:
            cursor.execute(statement, parameters)

        self.conn.commit()

    def get_locations(self):
        c = self.conn.cursor()
        c.execute("SELECT locId, locLon, locLat FROM tblLocations")
        locationss = c.fetchall()

        return locations


def build_location_database(location_db, location_file):
    from Utilities.shptools import shpReadShapeFile
    locations = []
    vertices, records = shpReadShapeFile(locationFile)
    for v, r in zip(vertices.values(), records):
        locLon, locLat = v[0]
        locId = str(r[0])
        locName = r[2]
        locCountry = r[4]
        locElev = r[7]
        locations.append((locId, locName, locLon, locLat,
                          locElev, locCountry, '', datetime.now()))

    return locations

def create_tables(conn):
    """
    Create a set of tables in the given database connection.

    :param conn: A valid sqlite3 database connection.

    """

    c = conn.cursor()
    try:
        c.execute(tblLocationsDef)
        c.execute(tblEventsDef)
        c.execute(tblWindSpeedDef)
        c.execute(tblHazardDef)
        c.execute(tblTracksDef)
    except sqlite3.OperationalError as e:
        log.warn("Warning: %s" % e.args[0])
        pass

    conn.commit()
    return conn

def setup(configFile, callback=None):
    """
    Set up the hazard database. This will create a hazard database in
    the output folders, create the tables for storing data and
    populate the tblLocations with details of all locations within the
    model domain.

    :param str configFile: Path to configuration file.

    """

    log.info("Generating hazard database")

    config = ConfigParser()
    config.read(configFile)

    outputPath = config.get('Output', 'Path')
    windfieldPath = pjoin(outputPath, 'windfield')
    trackPath = pjoin(outputPath, 'tracks')
    hazardPath = pjoin(outputPath, 'hazard')
    inputPath = './input'
    hazard_db = pjoin(outputPath, 'hazard.db')
    location_db = pjoin(inputPath, 'locations.db')

    if config.has_section('Region'):
        gridLimit = config.geteval('Region', 'gridLimit')

    log.info("Creating the hazard database")
    if os.path.isfile(hazard_db):
        log.info("Database file exists -- overwriting existing database")

    try:
        conn = sqlite3.connect(hazard_db)
    except sqlite3.Error as e:
        log.exception("Cannot create database: %s"% e.args[0])
        raise
    else:
        conn = create_tables(conn)
        conn.close()

    get_domain_locations(location_db, gridLimit, hazard_db)


def get_domain_locations(location_db, gridLimit, hazard_db):
    """
    Extract locations that lie within the model domain from the
    default locations database, and put them into a new database
    `hazard.db`.

    :param str location_db: Name of the location database file.
    :param dict gridLimit: The domain where the wind fields will be
                           generated. The :class:`dict` should contain
                           the keys :attr:`xMin`,:attr:`xMax`,
                           :attr:`yMin` and :attr:`yMax`. The *y*
                            variable bounds the latitude and the *x*
                            variable bounds the longitude.
    :param str hazard_db: Name of the output hazard database.



    """

    try:
        locdb = sqlite3.connect(location_db, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    except sqlite3.Error as e:
        log.exception("Cannot access the location database: %s" % e.args[0])
    else:
        c = locdb.cursor()
        c.execute(selectLocations, (gridLimit['xMin'], gridLimit['xMax'],
                                    gridLimit['yMin'], gridLimit['yMax']))

        locs = c.fetchall()
        locdb.close()

    try:
        hazdb = sqlite3.connect(hazard_db, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES )
    except sqlite3.Error as e:
        log.exception("Cannot access the hazard database: %s" % e.args[0])
    else:
        h = hazdb.cursor()
        h.executemany(insLocations, locs)
        hazdb.commit()
        hazdb.close()

    return

def windfield_attrs(ncfile):
    """
    Extract the required attributes from a netCDF file.

    :param str ncfile: Path to a valid netCDF file created by TCRM.

    :returns: A tuple containing the track filename, file modification date,
              TCRM version, minimum pressure and maximum wind, stored as
              global attributes in the netCDF file.
        
    """
    ncobj = Dataset(ncfile, 'r')
    trackfile = getattr(ncobj, 'track_file')
    trackfiledate = getattr(ncobj, 'track_file_date')
    trackfiledate = datetime.strptime(trackfiledate, '%Y-%m-%d %H:%M:%S')
    tcrm_version = getattr(ncobj, 'tcrm_version')
    trackfile = unicodedata.normalize("NFKD", trackfile).encode('utf-8', 'ignore')
    trackfile = os.path.basename(trackfile)
    tcrm_version = unicodedata.normalize("NFKD", tcrm_version).encode('utf-8', 'ignore')
    slpobj = ncobj.variables['slp']
    minslp = getattr(slpobj, 'actual_range')[0]

    vmaxobj = ncobj.variables['vmax']
    maxwind = getattr(vmaxobj, 'actual_range')[1]
    ncobj.close()
    return(trackfile, trackfiledate, tcrm_version, minslp, maxwind)
    
def generate_event_database(windfieldPath, hazard_db):
    """
    Populate tblEvents with details of the events. This stores the
    eventId, the wind field file, the corresponding track file, the
    maximum wind speed recorded in the event and the minimum sea level
    pressure.

    :param str windfieldPath: Path to the wind field files.
    :param str hazard_db: Path to the database file.

    """

    fileList = os.listdir(windfieldPath)
    files = [pjoin(windfieldPath, f) for f in fileList]
    files = [f for f in files if os.path.isfile(f)]

    try:
        hazdb = sqlite3.connect(hazard_db, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES )
    except sqlite3.Error as e:
        log.exception("Cannot access the hazard database: %s" % e.args[0])
        raise
    else:
        c = hazdb.cursor()

    params = []
    for n, f in enumerate(sorted(files)):
        log.debug("Processing {0}".format(f))
        si = os.stat(f)
        dtWindfieldFile = datetime.fromtimestamp(int(si.st_mtime))
        trackfile, dtTrackFile, tcrm_version, minslp, maxwind = windfield_attrs(f)
        params.append(("%06d"%n, os.path.basename(f), trackfile,
                       float(maxwind), float(minslp), dtTrackFile,
                       dtWindfieldFile, tcrm_version, "", datetime.now()))

    try:
        c.executemany(insEvents, params)
    except sqlite3.Error as e:
        log.exception("Cannot insert the records: %s" % e.args[0])
        raise
    else:
        hazdb.commit()
        hazdb.close()
    return
        
        
def process_events(windfieldPath, hazard_db):
    """
    Process the events (wind fields) for each location within the model domain. This will store the modelled wind speed (or the missing value) at each grid point, from each synthetic event.

    :param str windfieldPath: Path to the synthetic wind fields.
    :param str hazard_db: Name of the hazard database file.
    """

    fileList = os.listdir(windfieldPath)
    files = [pjoin(windfieldPath, f) for f in fileList]
    files = [f for f in files if os.path.isfile(f)]

    try:
        hazdb = sqlite3.connect(hazard_db, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES )
    except sqlite3.Error as e:
        log.exception("Cannot access the hazard database: %s" % e.args[0])
        raise

    c = hazdb.cursor()
    c.execute("SELECT locId, locLon, locLat FROM tblLocations")
    locs = c.fetchall()
    
    
    for n, f in enumerate(sorted(files)):
        log.info("Processing {0}".format(f))
        eventId = "%06d" % n
        ncobj = Dataset(f)
        lon = ncobj.variables['lon'][:]
        lat = ncobj.variables['lat'][:]
        
        vmax = ncobj.variables['vmax'][:]
        ua = ncobj.variables['ua'][:]
        va = ncobj.variables['va'][:]
        pmin = ncobj.variables['slp'][:]
        params = []
        
        for loc in locs:
            locId, locLon, locLat = loc
            i = find_index(lon, locLon)
            j = find_index(lat, locLat)
            locVm = vmax[j, i]
            locUa = ua[j, i]
            locVa = va[j, i]
            locPr = pmin[j, i]
            locParams = (locId, eventId, float(locVm), float(locUa),
                         float(locVa), float(locPr), " ", datetime.now())

            params.append(locParams)
        try:
            c.executemany(insWindSpeed, params)
        except sqlite3.Error as e:
            log.exception("Cannot access the hazard database: %s" % e.args[0])
            raise
        else:
            hazdb.commit()
            
    hazdb.close()

def process_hazard(hazardFile, hazard_db):
    """
    Update the hazard database with the return period information.

    :param str hazardFile: NetCDF file containing hazard (return period) data generated by TCRM.
    :param str hazard_db: Name of the hazard database file.

    """
    try:
        hazdb = sqlite3.connect(hazard_db, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES )
    except sqlite3.Error as e:
        log.exception("Cannot access the hazard database: %s" % e.args[0])
        raise

    c = hazdb.cursor()
    c.execute("SELECT locId, locLon, locLat FROM tblLocations")
    locs = c.fetchall()

    ncobj = Dataset(hazardFile)
    
    try:
        tcrm_version = getattr(ncobj, 'tcrm_version')
    except AttributeError:
        log.info("Missing tcrm_version attribute from {0}".format(hazardFile))
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

    params = []
    for k, year in enumerate(years):
        for loc in locs:
            locId, locLon, locLat = loc
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
        c.executemany(insHazard, params)
    except sqlite3.Error as e:
        log.exception("Cannot access the hazard database: %s" % e.args[0])
        raise
    else:
        hazdb.commit()

    hazdb.close()
