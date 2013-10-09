import sys
import numpy as np

from datetime import datetime
from matplotlib.dates import date2num, num2date
from scipy.interpolate import interp1d, splev, splrep

from Utilities.maputils import latLon2Azi
from Utilities.loadData import loadTrackFile

def interpolate(number, year, month, day, hour, minute, lon, lat, pressure, speed,
                bearing, windspeed, rmax, penv, delta, interpolation_type=None):
    """
    Interpolate the records in time to have a uniform time difference between
    records. Each of the input arrays represent the values for a single TC
    event. 
    
    :param number: `int` unique number for the TC to be interpolated.
    :param year: :class:`numpy.array` containing year of record
    :param month: :class:`numpy.array` containing month of record 
    :param day: :class:`numpy.array` containing day of record 
    :param hour: :class:`numpy.array` containing hour of record
    :param minute: :class:`numpy.array` containing minute of record
    :param lon: :class:`numpy.array` containing longitude of record
    :param lat: :class:`numpy.array` containing latitude of record
    :param pressure: :class:`numpy.array` containing central pressure of record
    :param speed: :class:`numpy.array` containing forward speed of TC of record
    :param bearing: :class:`numpy.array` containing bearing of record
    :param windspeed: :class:`numpy.array` containing maximum wind speed of
                      record
    :param rmax: :class:`numpy.array` containing radius of max winds of record
    :param penv: :class:`numpy.array` containing environmental pressure of
                 record
    :param delta: `float` time difference to interpolate the dataset to. Must be
                  positive.
    :param interpolation_type: Optional ['linear', 'akima'], specify the type
                               of interpolation used for the locations (i.e.
                               longitude and latitude) of the records.

    """
    
    day_ = [datetime(year[i], month[i], day[i], hour[i], minute[i])
            for i in xrange(year.size)]
    time_ = date2num(day_)
    dt_ = 24.0*np.diff(time_)
    dt = np.empty(hour.size, 'f')
    dt[1:] = dt_

    # Convert all times to a time after initial observation:
    timestep = 24.0*(time_ - time_[0])

    newtime = np.arange(timestep[0], timestep[-1]+.01, delta)
    newtime[-1] = timestep[-1]
    _newtime = (newtime/24.) + time_[0]
    newdates = num2date(_newtime)
    nid = number * np.ones(newtime.size)
    
    # FIXME: Need to address the issue when the time between obs is less 
    # than delta (e.g. two obs 5 hrs apart, but delta = 6 hrs). 

    if len(year) <= 2:
        # Use linear interpolation only (only a start and end point given):
        nLon = interp1d(timestep, lon, kind='linear')(newtime)
        nLat = interp1d(timestep, lat, kind='linear')(newtime)
        npCentre = interp1d(timestep, pressure, kind='linear')(newtime)
        npEnv = interp1d(timestep, penv, kind='linear')(newtime)
        nrMax = interp1d(timestep, rmax, kind='linear')(newtime)
        nwSpd = interp1d(timestep, windspeed, kind='linear')(newtime)
    else:
        if interpolation_type=='akima':
            # Use the Akima interpolation method:
            try:
                import _akima
            except ImportError:
                logger.exception( ("Akima interpolation module unavailable "
                                    " - default to scipy.interpolate") )
                nLon = splev(newtime, splrep(timestep, lon, s=0), der=0)
                nLat = splev(newtime, splrep(timestep, lat, s=0), der=0)
            else:
                nLon = _akima.interpolate(timestep, lon, newtime)
                nLat = _akima.interpolate(timestep, lat, newtime)
                
        elif interpolation_type=='linear':
            nLon = interp1d(timestep, lon, kind='linear')(newtime)
            nLat = interp1d(timestep, lat, kind='linear')(newtime)
        else:
            nLon = splev(newtime, splrep(timestep, lon, s=0), der=0)
            nLat = splev(newtime, splrep(timestep, lat, s=0), der=0)

        npCentre = interp1d(timestep, pressure, kind='linear')(newtime)
        nwSpd = interp1d(timestep, windspeed, kind='linear')(newtime)
        npEnv = interp1d(timestep, penv, kind='linear')(newtime)
        nrMax = interp1d(timestep, rmax, kind='linear')(newtime)



    bear_, dist_ = latLon2Azi(nLat, nLon, 1, azimuth=0)
    nthetaFm = np.zeros(newtime.size, 'f')
    nthetaFm[:-1] = bear_
    nthetaFm[-1] = bear_[-1]
    dist = np.zeros(newtime.size, 'f')
    dist[:-1] = dist_
    dist[-1] = dist_[-1]
    nvFm = dist/delta

    nYear = [date.year for date in newdates]
    nMonth = [date.month for date in newdates]
    nDay = [date.day for date in newdates]
    nHour = [date.hour for date in newdates]
    nMin = [date.minute for date in newdates]
    np.putmask(npCentre, npCentre > 10e+6, 0)
    
    return (nid, nYear, nMonth, nDay, nHour, nMin, nLon, nLat,
            npCentre, nvFm, nthetaFm, nwSpd, nrMax, npEnv)


def parseTracks(configFile, trackFile, source, delta, outputFile=None):
    """
    Load a track dataset, then interpolate to some time delta (given in
    hours). Events with only a single record are not altered.

    :type  configFile: string
    :param configFile: Configuration file containing settings that
                       describe the data source.

    :type  trackFile: string
    :param trackFile: Path to the input data source.

    :type  source: string
    :param source: Name of the data source. `configFile` must have a
                   corresponding section which contains options that
                   describe the data format.

    :type  delta: float
    :param delta: Time difference to interpolate the dataset to. Must be
                  positive.

    :type  outputFile: string
    :param outputFile: Path to the destination of output, if it is to
                       be saved.
                       
    """

    if delta < 0.0:
        raise ValueError("Time step for interpolation must be positive")
    
    i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe = \
                loadTrackFile(configFile, trackFile, source)

    results = []
    idx = np.flatnonzero(i)
    for n in xrange(len(idx)):
        if n != (len(idx) - 1):
            j = range(idx[n], idx[n+1])
        else:
            j = range(idx[n], len(i)-1)
            
        if len(j) == 1:
            # Save record directly:
            track = (n, y[j], m[j], d[j], h[j], mn[j], lon[j], lat[j],
                     p[j], s[j], b[j], w[j], r[j], pe[j])
        else:
            #Do an interpolation:
            track = interpolate(n, y[j], m[j], d[j], h[j], mn[j], lon[j],
                        lat[j], p[j], s[j], b[j], w[j], r[j],
                        pe[j], delta, 'linear')
            
        # Then save:
        results.append(track)
        
    newTracks = np.hstack([np.vstack(r) for r in results]).T
    
    if outputFile:
        header = ( 'TCID,Year,Month,Day,Hour,Minute,Longitude,'
                   'Latitude,CentralPressure,Speed,Bearing,Windspeed,'
                   'rMax,EnvPressure\n' )
                 
        fmt = '%i,%i,%i,%i,%i,%i,%8.3f,%8.3f,%7.2f,%6.2f,%6.2f,%6.2f,%6.2f,%7.2f'

        with open(outputFile, 'w') as fid:
            fid.write('%' + header)
            if len(newTracks) > 0:
                np.savetxt(fid, newTracks, fmt=fmt)

    return newTracks.T
