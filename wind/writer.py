"""
Saving the wind-field time-series to file

"""
import logging
import netCDF4
import affine
import numpy as np


class WriteFoliationCallback(object):
    """
    Incremental (by time step) recording of wind fields to NetCDF


    >>> import tempfile, datetime, xarray
    >>> with tempfile.NamedTemporaryFile(delete=False) as f:
    ...     res = 0.5
    ...     lat, lon = np.r_[-10:-5+.1:res], np.r_[110:115+.1:res]
    ...     gridLimit = dict(xMin=min(lon), yMin=min(lat),
    ...                      yMax=max(lat), xMax=max(lon))
    ...     callback = WriteFoliationCallback(f.name, gridLimit, res)
    ...     callback(datetime.datetime(2018,3,17),
    ...              *[np.ones((len(lat),len(lon)))]*4, lon=lon, lat=lat)
    ...     callback(datetime.datetime(2018,3,18), *[np.full((3,3),
    ...              fill_value=100)]*4, lon=lon[2:5], lat=lat[2:5])
    ...     callback(datetime.datetime(2018,3,20,12,0), *[np.ones((3,6))]*4,
    ...              lon=lon[0:6], lat=lat[2:5])
    ...     callback.ds.close()
    ...     x = xarray.open_dataset(f.name)
    >>> len(x.time)
    3
    >>> x.pressure.sum(dim=['lat','lon']).values.tolist()
    [121.0, 900.0, 18.0]
    >>> (x.gust_speed.isel(time=2, lon=slice(0,6),
    ...                    lat=slice(2,5)).values == 1).all()
    True

    Uses low-level NetCDF bindings to support windowed progressive writes.
    """

    def __init__(self, filename, gridLimit, resolution, margin=0,
                 maxchunk=256, wraps=None):
        """
        :param str filename: netCDF file for saving data to
        :param dict gridLimit: simulation domain extent
        :param float resolution: grid resolution in decimal degrees
        :param float margin: spatial extent over which the wind field is
                             calculated in decimal degrees
        :param int maxchunk: chunking size (for use in netCDF variable
                             creation)
        :param func wraps: Optional callback function (e.g. for
                           timeseries extraction)

        """

        logging.debug(
            "Preparing to record windfield evolution to {}".format(filename))

        self.callback = wraps

        def series(start, stop, inc=resolution):
            return np.linspace(start, stop, int(round((stop-start)/inc)) + 1)

        lat = series(gridLimit['yMin'] - margin, gridLimit['yMax'] + margin)
        lon = series(gridLimit['xMin'] - margin, gridLimit['xMax'] + margin)

        self.affine = ~(affine.Affine.translation(
                            xoff=gridLimit['xMin']-margin,
                            yoff=gridLimit['yMin']-margin) *
                        affine.Affine.scale(resolution))

        self.ds = root = netCDF4.Dataset(filename, mode='w')
        root.description = "Simulated Windfield Timeseries"

        # declare shapes
        root.createDimension('time', size=None)
        root.createDimension('lat', size=len(lat))
        root.createDimension('lon', size=len(lon))

        # coordinate variables:
        self.time = root.createVariable('time',
                                        datatype='f8',
                                        dimensions=('time',))
        self.lat = root.createVariable('lat',
                                       datatype='f4',
                                       dimensions=('lat',))
        self.lon = root.createVariable('lon',
                                       datatype='f4',
                                       dimensions=('lon',))

        self.lat[:] = lat
        self.lon[:] = lon

        # data variables:
        etc = dict(datatype='f4',
                   dimensions=('time', 'lat', 'lon'),
                   chunksizes=(1, min(len(lat), maxchunk),
                               min(len(lon), maxchunk)),
                   zlib=True)
        self.speed = root.createVariable('gust_speed', fill_value=0, **etc)
        self.Ux = root.createVariable('velocity_east', fill_value=0, **etc)
        self.Uy = root.createVariable('velocity_north', fill_value=0, **etc)
        self.P = root.createVariable('pressure', fill_value=np.NaN, **etc)

    def __call__(self, time, gust, Ux, Uy, P, lon, lat):
        """Save wind field layer for latest time step"""
        if self.callback:
            self.callback(time, gust, Ux, Uy, P, lon, lat)

        t = len(self.time)

        if not t:
            self.time.units = "days since " + time.strftime("%Y-%m-%d %H:%M")

        # convert window extent to slice indices
        origin = np.rint(self.affine * (lon[0], lat[0])).astype(int)
        opposite = np.rint(self.affine * (lon[-1], lat[-1])).astype(int)
        j, i = [slice(a, b+1) for a, b in zip(origin, opposite)]

        # append data
        self.time[t] = netCDF4.date2num(time, units=self.time.units)
        self.speed[t, i, j] = gust
        self.Ux[t, i, j] = Ux
        self.Uy[t, i, j] = Uy
        self.P[t, i, j] = P

        # save (flush) layer to file
        self.ds.sync()
