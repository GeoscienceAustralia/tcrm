"""
:mod:`definitions` -- table and statement definitions
=====================================================

.. module:: definitions
    :synopsis: Table definitions, insert statements and
               query statements for the database module.
.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

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
INSWINDSPEED = ("INSERT INTO tblWindSpeed "
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
