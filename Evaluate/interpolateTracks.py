import sys
import logging
import numpy as np

from datetime import datetime, timedelta
from matplotlib.dates import num2date
from scipy.interpolate import interp1d, splev, splrep

from Utilities.maputils import latLon2Azi
from Utilities.loadData import loadTrackFile, maxWindSpeed
from Utilities.track import Track, ncSaveTracks

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

TRACKFILE_COLS = ('Indicator', 'CycloneNumber', 'Year', 'Month',
                  'Day', 'Hour', 'Minute', 'TimeElapsed', 'Datetime',
                  'Longitude', 'Latitude', 'Speed', 'Bearing',
                  'CentralPressure', 'WindSpeed', 'rMax',
                  'EnvPressure')

TRACKFILE_FMTS = ('i', 'i', 'i', 'i',
                  'i', 'i', 'i', 'f', datetime,
                  'f', 'f', 'f', 'f', 'f',
                  'f', 'f', 'f', 'f')

TRACKFILE_OUTFMT = ('%i,%i,%i,%i,'
                    '%i,%i,%i,%8.3f,%s,%5.1f,'
                    '%8.3f,%8.3f,%6.2f,%6.2f,%7.2f,'
                    '%6.2f,%6.2f,%7.2f')

OUTPUT_COLS = ('CycloneNumber', 'Datetime', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'EnvPressure', 'rMax')

OUTPUT_FMTS = '%i,%s,%7.3f,%8.3f,%8.3f,%6.2f,%6.2f,%7.2f,%7.2f,%6.2f'

class Track2(object):
    def __init__(self, data):
        self.data = data
        self.trackId = None
        self.trackfile = None

    def __getattr__(self, key):
        if key.startswith('__') and key.endswith('__'):
            return super(Track, self).__getattr__(key)
        return self.data[key]


def interpolate(track, delta, interpolation_type=None):
    """
    Interpolate the records in time to have a uniform time difference between
    records. Each of the input arrays represent the values for a single TC
    event.

    :param track: :class:`Track` object containing all data for the track.
    :param delta: `float` time difference to interpolate the dataset to. Must be
                  positive.
    :param interpolation_type: Optional ['linear', 'akima'], specify the type
                               of interpolation used for the locations (i.e.
                               longitude and latitude) of the records.

    # FIXME: Need to address masking values - scipy.interpolate.interp1d
    handles numpy.ma masked arrays.
    """
    LOG.debug("Performing interpolation of TC track")
    if not hasattr(track, 'Datetime'):
        day_ = [datetime(*x) for x in zip(track.Year, track.Month,
                                          track.Day, track.Hour,
                                          track.Minute)]
    else:
        day_ = track.Datetime
    
    timestep = timedelta(delta/24.)
    try:
        time_ = np.array([d.toordinal() + (d.hour + d.minute/60.)/24.0
                      for d in day_], dtype=float)
    except AttributeError:
        import cftime
        if isinstance(day_[0], cftime.DatetimeJulian):
            day__ = [d._to_real_datetime() for d in day_]
            time_ = np.array([d.toordinal() + (d.hour + d.minute/60.)/24.
                              for d in day__], dtype=float)
        else:
            raise
    dt_ = 24.0 * np.diff(time_)
    dt = np.zeros(len(track.data), dtype=float)
    dt[1:] = dt_

    # Convert all times to a time after initial observation:
    timestep = 24.0*(time_ - time_[0])

    newtime = np.arange(timestep[0], timestep[-1] + .01, delta)
    newtime[-1] = timestep[-1]
    _newtime = (newtime / 24.) + time_[0]
    newdates = num2date(_newtime)
    newdates = np.array([n.replace(tzinfo=None) for n in newdates])

    if not hasattr(track, 'Speed'):
        idx = np.zeros(len(track.data))
        idx[0] = 1
        # TODO: Possibly could change `np.mean(dt)` to `dt`?
        track.WindSpeed = maxWindSpeed(idx, np.mean(dt), track.Longitude,
                                       track.Latitude, track.CentralPressure, 
                                       track.EnvPressure)
    # Find the indices of valid pressure observations:
    validIdx = np.where(track.CentralPressure < sys.maxsize)[0]

    # FIXME: Need to address the issue when the time between obs is less
    # than delta (e.g. only two obs 5 hrs apart, but delta = 6 hrs).

    if len(track.data) <= 3:
        # Use linear interpolation only (only a start and end point given):
        nLon = interp1d(timestep, track.Longitude, kind='linear')(newtime)
        nLat = interp1d(timestep, track.Latitude, kind='linear')(newtime)

        if len(validIdx) >= 2:
            npCentre = interp1d(timestep,
                                track.CentralPressure,
                                kind='linear')(newtime)
            nwSpd = interp1d(timestep,
                             track.WindSpeed,
                             kind='linear')(newtime)

        elif len(validIdx) == 1:
            # If one valid observation, assume no change and
            # apply value to all times
            npCentre = np.ones(len(newtime)) * track.CentralPressure[validIdx]
            nwSpd = np.ones(len(newtime)) * track.WindSpeed[validIdx]

        else:
            npCentre = np.zeros(len(newtime))
            nwSpd = np.zeros(len(newtime))

        npEnv = interp1d(timestep, track.EnvPressure, kind='linear')(newtime)
        nrMax = interp1d(timestep, track.rMax, kind='linear')(newtime)

    else:
        if interpolation_type == 'akima':
            # Use the Akima interpolation method:
            try:
                import akima
            except ImportError:
                LOG.exception(("Akima interpolation module unavailable "
                               " - default to scipy.interpolate"))
                nLon = splev(newtime, splrep(timestep, track.Longitude, s=0),
                             der=0)
                nLat = splev(newtime, splrep(timestep, track.Latitude, s=0),
                             der=0)

            else:
                nLon = akima.interpolate(timestep, track.Longitude, newtime)
                nLat = akima.interpolate(timestep, track.Latitude, newtime)

        elif interpolation_type == 'linear':
            nLon = interp1d(timestep, track.Longitude, kind='linear')(newtime)
            nLat = interp1d(timestep, track.Latitude, kind='linear')(newtime)

        else:
            nLon = splev(newtime, splrep(timestep, track.Longitude, s=0), der=0)
            nLat = splev(newtime, splrep(timestep, track.Latitude, s=0), der=0)

        if len(validIdx) >= 2:
            # No valid data at the final new time,
            # would require extrapolation:
            firsttime = np.where(newtime >= timestep[validIdx[0]])[0][0]
            lasttime = np.where(newtime <= timestep[validIdx[-1]])[0][-1]

            if firsttime == lasttime:
                # only one valid observation:
                npCentre = np.zeros(len(newtime))
                nwSpd = np.zeros(len(newtime))
                npCentre[firsttime] = track.CentralPressure[validIdx[0]]
                nwSpd[firsttime] = track.WindSpeed[validIdx[0]]

            else:
                npCentre = np.zeros(len(newtime))
                nwSpd = np.zeros(len(newtime))
                _npCentre = interp1d(timestep[validIdx],
                                     track.CentralPressure[validIdx],
                                     kind='linear')(newtime[firsttime:lasttime])

                _nwSpd = interp1d(timestep[validIdx],
                                  track.Speed[validIdx],
                                  kind='linear')(newtime[firsttime:lasttime])

                npCentre[firsttime:lasttime] = _npCentre
                nwSpd[firsttime:lasttime] = _nwSpd
                npCentre[lasttime] = _npCentre[-1]
                nwSpd[lasttime] = _nwSpd[-1]

        elif len(validIdx) == 1:
            npCentre = np.ones(len(newtime)) * track.CentralPressure[validIdx]
            nwSpd = np.ones(len(newtime)) * track.WindSpeed[validIdx]
        else:
            npCentre = np.zeros(len(newtime))
            nwSpd = np.zeros(len(newtime))

        npEnv = interp1d(timestep, track.EnvPressure, kind='linear')(newtime)
        nrMax = interp1d(timestep, track.rMax, kind='linear')(newtime)

    if len(nLat) >= 2:
        bear_, dist_ = latLon2Azi(nLat, nLon, 1, azimuth=0)
        nthetaFm = np.zeros(newtime.size, dtype=float)
        nthetaFm[:-1] = bear_
        nthetaFm[-1] = bear_[-1]
        dist = np.zeros(newtime.size, dtype=float)
        dist[:-1] = dist_
        dist[-1] = dist_[-1]
        nvFm = dist / delta

    else:
        nvFm = track.Speed[-1]
        nthetaFm = track.Bearing[-1]

    nYear = [date.year for date in newdates]
    nMonth = [date.month for date in newdates]
    nDay = [date.day for date in newdates]
    nHour = [date.hour for date in newdates]
    nMin = [date.minute for date in newdates]
    np.putmask(npCentre, npCentre > 10e+6, sys.maxsize)
    np.putmask(npCentre, npCentre < 700, sys.maxsize)

    newindex = np.zeros(len(newtime))
    newindex[0] = 1
    newTCID = np.ones(len(newtime)) * track.trackId[0]

    newdata = np.empty(len(newtime),
                        dtype={
                               'names': TRACKFILE_COLS,
                               'formats': TRACKFILE_FMTS
                               } )

    for key, val in zip(TRACKFILE_COLS,
                        [newindex, newTCID, nYear, nMonth, nDay, nHour, nMin,
                           newtime, newdates, nLon, nLat, nvFm, nthetaFm,
                           npCentre, nwSpd, nrMax, npEnv]):
        newdata[key] = val
    newtrack = Track(newdata)
    newtrack.trackId = track.trackId
    newtrack.trackfile = track.trackfile

    return newtrack

def saveTracks(tracks, outputFile):
    """
    Save the data to a TCRM-format track file (.csv format).

    :param tracks: A collection of `Track` objects to save to file.
    :param str outputFile: Path to destination file.

    """
    LOG.debug("Saving interpolated tracks to {0}".format(outputFile))
    output = []
    for track in tracks:
        r = [getattr(track, col).T for col in OUTPUT_COLS]
        output.append(r)
    data = np.hstack([r for r in output]).T

    with open(outputFile, 'w') as fp:
        fp.write('%' + ','.join(OUTPUT_COLS) + '\n')
        if len(data) > 0:
            np.savetxt(fp, data, fmt=OUTPUT_FMTS)

def parseTracks(configFile, trackFile, source, delta, outputFile=None,
                interpolation_type=None):
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

    :type  results: `list` of :class:`Track` objects containing the
                    interpolated track data

    """
    LOG.info("Interpolating tracks from {0}".format(trackFile))
    if delta < 0.0:
        raise ValueError("Time step for interpolation must be positive")

    if trackFile.endswith("nc"):
        from Utilities.track import ncReadTrackData
        tracks = ncReadTrackData(trackFile)
    else:
        tracks = loadTrackFile(configFile, trackFile, source)

    results = []

    for track in tracks:
        if len(track.data) == 1:
            results.append(track)
        else:
            newtrack = interpolate(track, delta, interpolation_type)
            results.append(newtrack)

    if outputFile:
        # Save data to file:
        ncSaveTracks(outputFile, results)


    return results
