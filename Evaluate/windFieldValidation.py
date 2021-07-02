import itertools
#import matplotlib.pyplot as plt
#from matplotlib import cm as cmap
import numpy as np

import xarray as xr
#import seaborn as sns

from wind import windmodels
from Utilities import metutils
from Utilities.maputils import bearing2theta, makeGrid, meshLatLon
from Utilities.parallel import attemptParallel


#sns.set_style('ticks', {'image.cmap':'coolwarm'})
#sns.set_context('poster')
#palette = [(1, 1, 1), (0.000, 0.627, 0.235), (0.412, 0.627, 0.235), (0.663, 0.780, 0.282),
#        (0.957, 0.812, 0.000), (0.925, 0.643, 0.016), (0.835, 0.314, 0.118),
#        (0.780, 0.086, 0.118)]
#cmap = sns.blend_palette(palette, as_cmap=True)

def polarGridAroundEye(lon, lat, margin=2, resolution=0.02):
    R, theta = makeGrid(lon, lat, margin, resolution)
    return R, theta

def meshGrid(lon, lat, margin=2, resolution=0.02):
    xgrid, ygrid = meshLatLon(lon, lat, margin, resolution)
    return xgrid, ygrid

def calculateWindField(lon, lat, pEnv, pCentre, rMax, vFm, thetaFm, beta,
                       profileType='powell', windFieldType='kepert'):

    pCentre = metutils.convert(pCentre, 'hPa', 'Pa')
    pEnv = metutils.convert(pEnv, 'hPa', 'Pa')
    vFm = metutils.convert(vFm, 'kmh', 'mps')
    thetaFm = bearing2theta(np.pi * thetaFm / 180.)
    thetaMax = 70.
    rmax = metutils.convert(rMax, 'km', 'm')
    cls = windmodels.profile(profileType)
    if profileType=="holland":
        profile = cls(lat, lon, pEnv, pCentre, rmax, beta)
    else:
        profile = cls(lat, lon, pEnv, pCentre, rmax)
    R, theta = polarGridAroundEye(lon, lat, 5.)
    gradV = profile.velocity(R*1000)
    cls = windmodels.field(windFieldType)
    windfield = cls(profile)
    Ux, Vy = windfield.field(R*1000, theta, vFm, thetaFm, thetaMax)

    surfV = np.sqrt(Ux*Ux+Vy*Vy)*1.268 # Gust conversion factor
    return gradV, surfV

"""
lat = np.arange(-30, -4, 2, dtype=float)
pc = np.arange(900, 991, 5, dtype=float)
pe = np.arange(995, 1016, dtype=float)
rm = np.arange(10, 91, 5, dtype=float)
vfm = np.arange(0, 51, 5, dtype=float)
gwind = np.zeros((len(lat), len(pc), len(pe), len(rm), len(vfm)))
swind = np.zeros((len(lat), len(pc), len(pe), len(rm), len(vfm)))
it = np.nditer(gwind, flags=['multi_index'])
nn = gwind.size
print(nn)

lon = 120.
thetaFm = 70
beta = 1.6
profileType = "powell"
blmodel = "kepert"
i = 0



for x in it:
    il, ic, ip, ir, iv = it.multi_index
    gradV, surfV = calculateWindField(lon, lat[il], pe[ip], pc[ic],
                            rm[ir], vfm[iv], thetaFm, beta,
                            profileType=profileType,
                            windFieldType=blmodel)
    gwind[it.multi_index] = np.max(gradV)
    swind[it.multi_index] = np.max(surfV)
    i += 1
    print(f"{100*i/nn:0.4f} %")

coords = [
    ("latitude", lat, dict(long_name="Latitude",
                           units="degrees_south")),
    ("pcentre", pc, dict(long_name="Central pressure",
                         units="hPa")),
    ("penv", pe, dict(long_name="Environmental pressure",
                      units="hPa")),
    ("rmax", rm, dict(long_name="Radius to maximum winds",
                      units="km")),
    ("vfm", vfm, dict(long_name="Forward speed",
                      units="km/h"))
]

dims = ["latitude", 'pcentre', 'penv', 'rmax', 'vfm']
gattrs = {
    "long_name": "Gradient level wind speed",
    "profile": profileType,
    "blmodel": blmodel,
    "description": "maximum gradient level wind speed",
    "units": "m s-1",
    }
sattrs = {
    "long_name": "Surface wind speed",
    "profile": profileType,
    "blmodel": blmodel,
    "description": "maximum 0.2-s wind gust",
    "units": "m s-1",
    }


gda = xr.DataArray(gwind, dims=dims, coords=coords, attrs=gattrs)
sda = xr.DataArray(swind, dims=dims, coords=coords, attrs=sattrs)
ds = xr.Dataset()
ds['gradwind'] = gda
ds['surfwind'] = sda
ds.to_netcdf("output.nc")
"""

def balanced(iterable):
    """
    Balance an iterator across processors.

    This partitions the work evenly across processors. However, it
    requires the iterator to have been generated on all processors
    before hand. This is only some magical slicing of the iterator,
    i.e., a poor man version of scattering.
    """
    P, p = MPI.COMM_WORLD.size, MPI.COMM_WORLD.rank
    return itertools.islice(iterable, p, None, P)

def run():
    lat = np.arange(-30, -4, 2, dtype=float)
    pc = np.arange(900, 991, 5, dtype=float)
    pe = np.arange(995, 1016, dtype=float)
    rm = np.arange(10, 91, 5, dtype=float)
    vfm = np.arange(0, 51, 5, dtype=float)
    gwind = np.zeros((len(lat), len(pc), len(pe), len(rm), len(vfm)))
    swind = np.zeros((len(lat), len(pc), len(pe), len(rm), len(vfm)))
    it = np.nditer(gwind, flags=['multi_index'])
    nn = gwind.size
    #print(nn)

    lon = 120.
    thetaFm = 70
    beta = 1.6
    profileType = "powell"
    blmodel = "kepert"
    i = 0

    # Attempt to start the track generator in parallel
    global MPI
    MPI = attemptParallel()
    comm = MPI.COMM_WORLD

    status = MPI.Status()
    worktag = 0
    resulttag = 1
    idx = [it.multi_index for x in it]
    
    if (comm.rank == 0) and (comm.size > 1):
        w = 0
        p = comm.size -1
        for d in range(1, comm.size):
            print(w)
            if w < len(idx):
                comm.send(idx[w], dest=d, tag=worktag)
                w += 1
            else:
                comm.send(None, dest=d, tag=worktag)
                p = w

        terminated = 0

        while terminated < p:
            try:
                result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            except Exception:
                pass

            d = status.source
            if result:
                gV, sV, workidx = result
                gwind[workidx] = gV
                swind[workidx] = sV
                #gwind[idx[w]], swind[idx[w]] = result

            if w < len(idx):
                comm.send(idx[w], dest=d, tag=worktag)
                w += 1
            else:
                comm.send(None, dest=d, tag=worktag)
                terminated += 1

    elif (comm.rank != 0) and (comm.size > 1):
        while True:
            workidx = comm.recv(source=0, tag=worktag, status=status)
            if workidx is None:
                break
            il, ic, ip, ir, iv = workidx
            print(f"Processing {workidx}")
            gradV, surfV = calculateWindField(lon, lat[il], pe[ip], pc[ic],
                                rm[ir], vfm[iv], thetaFm, beta,
                                profileType=profileType,
                                windFieldType=blmodel)
            results = (np.max(np.abs(gradV)), np.max(surfV), workidx)
            comm.send(results, dest=0, tag=resulttag)

    elif (comm.rank == 0) and (comm.size == 1):
        for x in idx:
            il, ic, ip, ir, iv = x
            print(lat[il], pc[ic], pe[ip], rm[ir], vfm[iv])
            gradV, surfV = calculateWindField(lon, lat[il], pe[ip], pc[ic],
                                              rm[ir], vfm[iv], thetaFm, beta,
                                              profileType=profileType,
                                              windFieldType=blmodel)
            gwind[x] = np.max(np.abs(gradV))
            swind[x] = np.max(surfV)

    comm.barrier()

    coords = [
        ("latitude", lat, dict(long_name="Latitude",
                            units="degrees_south")),
        ("pcentre", pc, dict(long_name="Central pressure",
                            units="hPa")),
        ("penv", pe, dict(long_name="Environmental pressure",
                        units="hPa")),
        ("rmax", rm, dict(long_name="Radius to maximum winds",
                        units="km")),
        ("vfm", vfm, dict(long_name="Forward speed",
                        units="km/h"))
    ]

    dims = ["latitude", 'pcentre', 'penv', 'rmax', 'vfm']
    gattrs = {
        "long_name": "Gradient level wind speed",
        "profile": profileType,
        "blmodel": blmodel,
        "description": "maximum gradient level wind speed",
        "units": "m s-1",
        }
    sattrs = {
        "long_name": "Surface wind speed",
        "profile": profileType,
        "blmodel": blmodel,
        "description": "maximum 0.2-s wind gust",
        "units": "m s-1",
        }

    if comm.rank == 0:
        gda = xr.DataArray(gwind, dims=dims, coords=coords, attrs=gattrs)
        sda = xr.DataArray(swind, dims=dims, coords=coords, attrs=sattrs)
        ds = xr.Dataset()
        ds['gradwind'] = gda
        ds['surfwind'] = sda
        ds.to_netcdf("output.nc")

    MPI.Finalize()

if __name__ == '__main__':
    print("Starting")
    global MPI, comm
    print("Initialiszing MPI")
    MPI = attemptParallel()
    #import atexit
    #atexit.register(MPI.Finalize)
    comm = MPI.COMM_WORLD

    print("Executing run()")
    run()

    #MPI.Finalize()
