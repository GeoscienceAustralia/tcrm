"""
:mod:`metutils` -- perform basic meteorological calculations
============================================================

.. module:: metutils
    :synopsis: perform basic meteorological calculations.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import math
import numpy as np
import numpy.ma as ma

#Define constants
gPressureUnits = "hPa"
gApproxPressure = 101.325
gEps = 0.622

def elevToAirPr(elev, units_ap=gPressureUnits):
    """
    Approximate air pressure in hectopascals (mb) at a given elevation.

    :param float elev: Elevation above sea level (metres).
    :param units_ap: Units of air pressure (default ``"hPa"``).
    :type units_ap: str

    :returns: Approximate air pressure at the given elevation
              in the specified or default units.
    :rtype: float

    """

    ap = gApproxPressure
    if elev > 0:
        ap = gApproxPressure * np.exp(-0.0001184 * elev)

    ap = convert(ap, 'kPa', units_ap)
    return ap

def vapPrToDewPoint(vp, units_vp=gPressureUnits):
    """
    Calculate dew point from vapour pressure (in kPa)

    :param float vp: Input vapour pressure.
    :param str units_vp: Input units (default ``gPressureUnits``)

    :returns: dew point temperature (in degrees Kelvin)
    :rtype: float

    """
    vp = convert(vp, units_vp, "kPa")
    t_dp = (116.9 + (237.3 * np.log(vp))) / (16.78 - np.log(vp))
    return t_dp

def dewPointToVapPr(t_dp, units_vp=gPressureUnits):
    """
    Calculate vapour pressure from dew point temperature.

    :param float t_dp: Dew point temperature (degrees Celsius)
    :param str units_vp: Output units (default ``gPressureUnits``).

    :returns: Vapour pressure of water content.
    :rtype: float

    """
    vp = np.exp((16.78 * t_dp - 116.9)/(t_dp + 237.3))
    vp = convert(vp, 'kPa', units_vp)
    return vp

def wetBulbGlobeTemp(t_dp, temp):
    """
    Calculate Wet Bulb Globe Temperature from Dew Point Temperature
    and Dry Bulb Temperature (same as air temperature).

    :param float t_dp: Dew point temperature (degees Celsius).
    :param float temp: Air temperature (degrees Celsius).

    :returns: Wet bulb globe temperature (degrees Celsius).
    :rtype: float

    """
    vp = dewPointToVapPr(t_dp, 'kPa')
    wbgt = 0.567 * temp + 0.393 * vp + 3.94
    return wbgt

def wetBulbToDewPoint(db, wb, elev=0):
    """
    Calculate Dew Point from dry bulb and wet bulb temperatures.

    :param float db: Dry bulb temperature (degrees Celsius).
    :param float wb: Wet bulb temperature (degrees Celsius).
    :param float elev: Optional elevation of the observation (metres).

    :returns: Dew point temperature (degrees Kelvin).
    :rtype: float

    """
    # Calculate vapour pressure
    vp = wetBulbToVapPr(db, wb, elev, 'kPa')
    # Dew point
    t_dp = vapPrToDewPoint(vp, 'kPa')
    return t_dp

def wetBulbToVapPr(db, wb, elev, units_vp=gPressureUnits):
    """
    Calculate vapour pressure from dry bulb and wet bulb temperatures,
    and optional elevation.

    :param float db: Dry bulb temperature (degrees Celsius).
    :param float wb: Wet bulb temperature (degrees Celsius).
    :param float elev: Elevation (metres).
    :param str units_vp: Output units (default ``gPressureUnits``)

    :returns: Vapour pressure.
    :rtype: float

    """
    if (wb > db):
        # Reality check. Wet bulb can't be greater than dry bulb
        wb = db

    # Get saturation vapour pressure at wet bulb temperature, in kPa.
    sat_vp_wb = satVapPr(wb, 'kPa')

    # Conversion factor, $cfA
    cfA = 0.00066 * (1 + 0.00115 * wb)

    # Air pressure
    ap = elevToAirPr(elev, 'kPa')

    # Vapour pressure (partial pressure of water vapour in kPa)
    vp = sat_vp_wb - (cfA * ap * (db - wb))
    vp = convert(vp, 'kPa', units_vp)
    return vp

def satVapPr(temp, units_vp=gPressureUnits):
    """
    Saturation vapour pressure from temperature in degrees celsius.

    :param float temp: Temperature (degrees celsius).
    :param str units_vp: Units of the vapour pressure to return.
                         Default is ``gPressureUnits``.

    :returns: saturation vapour pressure in the specified or default units.

    .. note::
       Calculation is in kPa with a conversion at the end to the
       required units.


    Example::

        >>> from metutils import satVapPr
        >>> satVapPr(25.)
        31.697124349060619

    """
    vp = np.exp(((16.78 * temp) - 116.9) / (temp + 237.3))
    vp = convert(vp, 'kPa', units_vp)

    return vp

def vapPrToRH(vp, sat_vp):
    """
    Calculate relative humidity from vapour pressure and saturated
    vapour pressure.

    :param float vp: Vapour pressure (hPa)
    :param float sat_vp: Saturation vapour pressure (hPa)

    :returns: Relative humidity (%)

    Example::

        >>> from metutils import vapPrToRH
        >>> vapPrToRH(10., 30.)
        33.33333333333

    """
    if (sat_vp == 0):
        rh = 100
    else:
        rh = (vp / sat_vp) * 100.

    # Any out of bounds value is converted to a boundary value
    if rh > 100:
        rh = 100
    elif rh < 0:
        rh = 0
    return rh

def wetBulbToRH(t_db, t_wb, elev):
    """
    Calculate relative humidity from dry bulb and wet bulb temperatures,
    and optional elevation.

    :param float t_db: Dry bulb temperature (degrees Celsius).
    :param float t_wb: Wet bulb temperature (degrees Celsius).
    :param float elev: Elevation (metres).

    :returns: Relative humidity (%).

    Example::

        >>> from metutils import wetBulbToRH
        >>> wetBulbToRH(25., 10., 0)
        6.747686
        >>> wetBulbToRH(25., 20., 0)
        63.01954

    """

    # Calculate vapour pressure
    vp = wetBulbToVapPr(t_db, t_wb, elev)
    # Calculate the saturated vapour pressure at the dry bulb temperature
    sat_vp = satVapPr(t_db, 'kPa')

    # Calculate the relative humidity
    rh = vapPrToRH(vp, sat_vp)
    return rh

def dewPointToRH(t_dry, t_dew):
    """
    Calculate relative humidity from dry bulb and dew point (in degrees
    Celsius)

    :param float t_dry: Dry bulb temperature (degrees Celsius).
    :param float t_dew: Dew point temperature (degrees Celsius).

    :returns: Relative humidity (%)
    :rtype: float

    """
    vap_t_dry = 6.11 * (10 ** ((7.5 * t_dry) / (237.3 + t_dry)))
    vap_t_dew = 6.11 * (10 ** ((7.5 * t_dew) / (237.3 + t_dew)))
    rh = (vap_t_dew / vap_t_dry) * 100.
    # Any out of bounds value is converted to an undefined value
    if (rh > 100 or rh < 0):
        rh = None
    return rh

def rHToDewPoint(rh, t_dry):
    """
    Calculate dew point from relative humidity and dry bulb (in degrees
    Celsius).

    :param float rh: Relative humidity (%).
    :param float t_dry: Dry bulb temperature (degrees Celsius).

    :returns: Dew point temperture (degrees Celsius).
    :rtype: float
    """
    vap_t_dry = satVapPr(t_dry, 'kPa')
    vap_t_dew = (rh/100.0) * vap_t_dry
    if (vap_t_dew > 0):
        t_dew = vapPrToDewPoint(vap_t_dew, "kPa")

    # Any out of bounds value is converted to an undefined value
    if (t_dew > t_dry):
        t_dew = None
    return t_dew

def vapPrToMixRat(es, prs):
    """
    Calculate mixing ratio from vapour pressure
    In this function, we (mis)use the symbol es for vapour pressure,
    when it correctly represents saturation vapour pressure.
    The function can be used for both.

    :param float es: Vapour pressure.
    :param float prs: Air pressure.

    :returns: Mixing ratio.
    :rtype: float

    """
    rat = gEps * es / (prs - es)
    return rat

def mixRatToVapPr(rat, prs):
    """
    Calculate vapour pressure from mixing ratio.

    :param float rat: Mixing ratio.
    :param float prs: Air pressure.

    :returns: Vapour pressure.
    :rtype: float
    """
    es = rat * prs / (gEps + rat)
    return es

def vapPrToSpHum(es, prs):
    """
    Convert vapour pressure to specific humidity.

    :param float es: Vapour pressure.
    :param float prs: Air pressure.

    :returs: Specific humidity.
    :rtype: float
    """
    q = gEps * es / prs
    return q

def spHumToMixRat(q, units="gkg"):
    """
    Calculate mixing ratio from specific humidity.
    Assumes the input specific humidity variable is in units
    of g/kg.

    :param float q: Specific humidity (any units, assumed g/kg).
    :param str units: Units of specific humidity, default g/kg.

    :returns: Mixing ratio (kg/kg).
    :rtype: float

    """
    q = convert(q, units, "kgkg")

    rat = gEps * q / (gEps - q)
    return rat

def rHToMixRat(rh, tmp, prs, tmp_units="C"):
    """
    Calculate mixing ratio from relative humidity, temperature and pressure.

    :param float rh: Relative humidity (%).
    :param float tmp: Temperature (any units, default degrees Celsius).
    :param float prs: Air pressure (hPa).
    :param str tmp_units: Air temperature units (default degrees Celsius).

    :returns: Mixing ratio (g/kg).
    :rtype: float

    """
    es = satVapPr(convert(tmp, tmp_units, "C"))
    e = (rh / 100.) * es
    rat = vapPrToMixRat(e, prs)
    return rat

def spHumToRH(q, tmp, prs):
    """
    Calculate relative humidity from specific humidity, temperature and
    pressure.

    :param float q: Specific humidity (g/kg).
    :param float tmp: Temperature (degrees Celsius).
    :param float prs: Air pressure (hPa).

    :returns: Relative humidity (%)
    :rtype: float

    """
    es = satVapPr(tmp)
    qs = gEps * es / prs
    rh = 100. * q / qs
    return rh

def coriolis(lat):
    """
    Calculate the Coriolis factor (f) for a given latitude (degrees).

    :param lat: Latitude (degrees).
    :type  lat: :class:`numpy.ndarray` or scalar float

    :returns: Coriolis factor
    :rtype: :class:`numpy.ndarray` or scalar float

    """
    omega = 2 * math.pi / 86400.
    f = 2 * omega * np.sin(np.radians(lat))
    return f

def convert(value, inunits, outunits):
    """
    Convert value from input units to output units.

    :param value: Value to be converted
    :param str inunits: Input units.
    :param str outunits: Output units.

    :returns: Value converted to ``outunits`` units.

    """
    startValue = value
    value = ma.array(value, dtype=float)
    if inunits == outunits:
        # Do nothing:
        return value
    if inunits == 'kmh':
        inunits = 'kph'
    if inunits == "m/s":
        inunits = "mps"

    # Speeds:
    mps = {"kph":3.6, "kts":1.944, "mph":2.2369}
    mph = {"kph":1.60934, "kts":0.86898, "mps":0.44704}
    kph = {"kts":0.539957, "mps":0.2777778, "mph":0.621371}
    kts = {"kph":1.852, "mps":0.5144, "mph":1.15}

    # Temperatures:
    C = {"F":1.8, "K":1.}
    F = {"C":0.5556}
    K = {"C":1.}

    # Pressures:
    kPa = {"hPa":10., "Pa":1000., "inHg":0.295299831, "mmHg":7.500615613,
           "Pascals":1000.}
    hPa = {"kPa":0.1, "Pa":100., "inHg":0.02953, "mmHg":0.750061561,
           "Pascals":100.}
    Pa = {"kPa":0.001, "hPa":0.01, "inHg":0.0002953, "mmHg":0.007500616,
          "Pascals":1.0}
    inHg = {"kPa":3.386388667, "hPa":33.863886667, "Pa":3386.388666667,
            "mmHg":25.4}
    mmHg = {"kPa":0.13332239, "hPa":1.3332239, "Pa":133.32239, "inHg":0.0394}
    pascals = {"kPa":0.001, "hPa":0.01, "inHg":0.0002953, "mmHg":0.007500616,
          "Pa":1.0}

    # Lengths:
    km = {"m":1000., "mi":0.621371192, "deg":0.00899886, "nm":0.539957,
          "rad":0.0001570783}
    deg = {"km":111.1251, "m":111125.1, "mi":69.0499358, "nm":60.0,
           "rad":math.pi/180.}
    m = {"km":0.001, "mi":0.000621371, "deg":0.00000899886, "nm":0.000539957,
         "rad":0.0000001570783}
    mi = {"km":1.60934, "m":1609.34, "deg":0.014482}
    nm = {"km":1.852, "m":1852, "deg":0.01666, "rad":math.pi/10800.}
    rad = {"nm":10800./math.pi, "km":6366.248653, "deg":180./math.pi}

    # Mixing ratio:
    gkg = {"kgkg":0.001}
    kgkg = {"gkg":1000}

    convert = {"mps":mps,
               "mph":mph,
               "kph":kph,
               "kts":kts,
               "kPa":kPa,
               "hPa":hPa,
               "Pa":Pa,
               "Pascals":pascals,
               "inHg":inHg,
               "mmHg":mmHg,
               "C":C,
               "F":F,
               "K":K,
               "km":km,
               "m":m,
               "deg":deg,
               "mi":mi,
               "nm":nm,
               "rad":rad,
               "gkg":gkg,
               "kgkg":kgkg}

    # Additions required before multiplication:
    convert_pre = {"F":{"C":-32.}}
    # Additions required after multiplication:
    convert_post = {"C":{"K":273., "F":32.},
                    "K":{"C":-273.}}

    if inunits in convert_pre:
        if outunits in convert_pre[inunits]:
            value += convert_pre[inunits][outunits]

    if inunits in convert:
        if outunits in convert[inunits]:
            value = value * convert[inunits][outunits]

    if inunits in convert_post:
        if outunits in convert_post[inunits]:
            value += convert_post[inunits][outunits]

    return value

def vapour(temp):
    """
    Determine saturation vapour pressure, given temperature).

    :param float temp: Air temperature (degrees Celsius)

    :returns: Saturation vapour pressure (kPa).
    :rtype: float

    """
    vp = 6.112 * np.exp(17.67 * temp/(243.5 + temp))
    return vp

def genesisPotential(zeta, rh, vmax, shear):
    """
    Calculate genesis potential index
    """
    gpi = np.power(abs((10 ** 5) * zeta), 1.5) * \
          ((rh / 50.) ** 3) * \
          ((vmax/70.) ** 3) / \
          ((1. + 0.1 * shear) ** 2)
    return gpi

def dewPointToWetBulb(T, Td, pressure):
    """
    Calculate wet bulb temperature from dry bulb temperature,
    dew point temperature and air pressure.

    Based on the code underpinning the form at:
    http://www.srh.noaa.gov/epz/?n=wxcalc_dewpoint

    :param float T: Dry bulb temperature (degrees Celcius).
    :param float Td: Dew point temperature (degrees Celcius).
    :param float pressure: Air pressure (hPa).

    :returns: Wet bulb temperature (degrees Celcius).
    :rtype: float

    TODO: Add unit tests.
    """

    if Td > T:
        raise ValueError("Dew point cannot be higher than dry bulb temperature")

    Edifference = 1
    incr = 10.
    prevsign = 1
    Tw = Td
    Es = satVapPr(Td, 'hPa')
    while (abs(Edifference) > 0.05):
        Ewguess = satVapPr(Tw, 'hPa')
        Eguess = Ewguess - pressure * (T - Tw) * \
                 0.00066 * (1. + (0.00115 * Tw))
        Edifference = Es - Eguess

        if (Edifference == 0):
            break
        else:
            if (Edifference < 0.):
                cursign = -1
                if (cursign != prevsign):
                    prevsign = cursign
                    incr = incr/10.
            else:
                cursign = 1
                if (cursign != prevsign):
                    prevsign = cursign
                    incr = incr/10.

        if (abs(Edifference) <= 0.05):
            break
        else:
            Tw = Tw + incr * prevsign

    return np.round(Tw, decimals=3)

