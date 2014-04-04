"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Title: met.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-11-28
Description: Contains functions to perform met calculations

Members:
elevToAirPr(elev <, units_ap>) :
vapPrToDewPoint(vp) :
dewPointToVapPr(t_dp <, units_vp>) :
wetbulbgt(t_dp, temp) :
wetBulbToDewPoint(db, wb, elev) :
wetBulbToVapPr(db, wb, elev <, units_vp>) :
satVapPr(temp <, units_vp>) :
vapPrToRH(vp, sat_vp) :
wetBulbToRH(t_db, t_wb, elev) :
dewPointToRH(t_dry, t_dew) :
rHToDewPoint(rh, t_dry) :
vapPrToMixRat(es, prs) :
rHToMixRat(rh, tmp, prs):
coriolis(latitude) :
    Returns the coriolis parameter (commonly called 'f') for a given
    latitude convert(value, input, output):
    Converts 'value' from 'input' units to 'output' units. Contains
    conversions for length, speed, pressure and temperature.

Version: $Rev: 642 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-09-11
Modification: Added functions for converting atmospheric moisture variables

SeeAlso: (related programs)
Constraints:

$Id: metutils.py 642 2012-02-21 07:54:04Z nsummons $
"""

"""
Contents:
coriolis(lat) :
    Calculates Coriolis parameter (f) from latitude (in degress)

convert(value, input, output) :
    Convert value from input units to gPressureUnitsoutput units.
    'input' and 'output' are strings, 'value' is a float

"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
import math
import numpy
from numpy import radians
__version__ = '$Id: metutils.py 642 2012-02-21 07:54:04Z nsummons $'

#Define constants
gPressureUnits = "hPa"
gApproxPressure = 101.325
gEps = 0.622

def elevToAirPr(elev, units_ap=gPressureUnits):
    """elevToAirPr(elev <, units_ap>)
    Approximate air pressure in hectopascals (mb)
    Input: elevation above sea level (metres)
    Output: approx. air pressure in the specified or default units.
    Note: Calculation is in kPa with a conversion at the end to the
          required units.
    """
    ap = gApproxPressure
    if elev > 0:
        ap = gApproxPressure * numpy.exp(-0.0001184 * elev)

    ap = convert(ap, 'kPa', units_ap)
    return ap

def vapPrToDewPoint(vp, units_vp=gPressureUnits):
    """vapPrToDewPoint(vp)
    Calculate Dew Point from vapour pressure (in kPa)
    """
    vp = convert(vp, units_vp, "kPa")
    t_dp = (116.9 + (237.3*numpy.log(vp))) / (16.78 - numpy.log(vp))
    return t_dp

def dewPointToVapPr(t_dp, units_vp=gPressureUnits):
    """dewPointToVapPr(t_dp <, units_vp>)
    Calculate vapour pressure (in default units (hPa) or specified
    units) from dew point temperature
    """
    vp = numpy.exp((16.78*t_dp-116.9)/(t_dp+237.3))
    vp = convert(vp, 'kPa', units_vp)
    return vp

def wetBulbGlobeTemp(t_dp, temp):
    """wetBulbGlobeTemp(t_dp, temp)
    Calculate Wet Bulb Globe Temperature from Dew Point Temperature
    and Dry Bulb Temperature (same as air temperature).
    Returns null if dew point or temp not defined
    """
    vp = dewPointToVapPr(t_dp, 'kPa');
    wbgt = 0.567*temp + 0.393*vp + 3.94;
    return wbgt

def wetBulbToDewPoint(db, wb, elev=0):
    """wetBulbToDewPoint(db, wb, elev)
    Calculate Dew Point from dry bulb and wet bulb temperatures
    """
    # Calculate vapour pressure
    vp = wetBulbToVapPr(db, wb, elev, 'kPa');
    # Dew point
    t_dp = vapPrToDewPoint(vp, 'kPa')
    return t_dp

def wetBulbToVapPr(db, wb, elev, units_vp=gPressureUnits):
    """wetBulbToVapPr(db, wb, elev, units_vp=gPressureUnits )
    Calculate vapour pressure from dry bulb and wet bulb temperatures,
    and optional elevation.
    input: dry bulb and wet bulb temps in degrees centigrade, elevation
    in metres (elevation is optional)
    output: vapour pressure
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
    vp = sat_vp_wb - (cfA*ap*(db - wb))

    return vp

def satVapPr(temp, units_vp=gPressureUnits):
    """satVapPr(temp, units_vp=gPressureUnits)
    Saturation vapour pressure from temperature in degrees celsius.
    Input: temperature (degrees celsius)
    Output: saturation vapour pressure in the specified or default units.
    Note: Calculation is in kPa with a conversion at the end to the
          required units.
    """
    vp = numpy.exp(((16.78 * temp) - 116.9) / (temp+237.3))
    vp = convert(vp, 'kPa', units_vp)

    return vp

def vapPrToRH(vp, sat_vp):
    """vapPrToRH(vp, sat_vp)
    Calculate relative humidity from vapour pressure and saturated
    vapour pressure
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
    """wetBulbToRH(t_db, t_wb, elev)
    Calculate relative humidity from dry bulb and wet bulb temperatures,
    and optional elevation.
    input: dry bulb and wet bulb temps in degrees centigrade, elevation
           in metres (elevation is optional)
    output: relative humidity
    """
    # Calculate vapour pressure
    vp = wetBulbToVapPr(t_db, t_wb, elev)
    # Calculate the saturated vapour pressure at the dry bulb temperature
    sat_vp = satVapPr(t_db, 'kPa')

    # Calculate the relative humidity
    rh = vapPrToRH(vp, sat_vp)
    return rh

def dewPointToRH(t_dry, t_dew):
    """dewPointToRH(t_dry, t_dew)
    Calculate relative humidity from dry bulb and dew point (in degrees
    Celsius)

    Calculation take from Dave Williamson's Access code:
        Val([sTempDryBulb]) AS Tdry,
        Val([sDewPoint]) AS Tdew,
        6.11*10^((7.5*[Tdry])/(273.3+[Tdry])) AS VapT,
        6.11*10^((7.5*[Tdew])/(273.3+[Tdew])) AS VapTd,
        CInt(100*[VapTd]/[VapT]) AS RH
    """
    vap_t_dry = 6.11 * (10 ** ((7.5 * t_dry) / (237.3 + t_dry)))
    vap_t_dew = 6.11 * (10 ** ((7.5 * t_dew) / (237.3 + t_dew)))
    rh = (vap_t_dew / vap_t_dry)*100.
    # Any out of bounds value is converted to an undefined value
    if (rh > 100 or rh < 0):
        rh = None
    return rh

def rHToDewPoint(rh, t_dry):
    """rHToDewPoint(rh, t_dry)
    Calculate dew point from relative humidity and dry bulb (in degrees
    Celsius)
    """
    vap_t_dry = satVapPr(t_dry, 'kPa')  # 6.11 * (10 ** ((7.5 * t_dry ) / (237.3 + t_dry)))
    vap_t_dew = (rh/100.0) * vap_t_dry
    if (vap_t_dew > 0):
        t_dew = vapPrToDewPoint(vap_t_dew, "kPa")

    # Any out of bounds value is converted to an undefined value
    if (t_dew > t_dry):
        t_dew = None
    return t_dew

def vapPrToMixRat(es, prs):
    """vapPrToMixRat(es, prs)
    Calculate mixing ratio from vapour pressure
    In this function, we (mis)use the symbol es for vapour pressure,
    when it correctly represents saturation vapour pressure.
    The function can be used for both
    """
    rat = gEps*es/(prs-es)
    return rat

def mixRatToVapPr(rat, prs):
    """mixRatToVapPr(rat, prs):
    Calculate vapour pressure from mixing ratio
    """
    es = rat*prs/(gEps+rat)
    return es

def vapPrToSpHum(es, prs):
    """vapPrToSpHum(es, prs)
    Convert vapour pressure to specific humidity
    """
    q = gEps*es/prs
    return q

def spHumToMixRat(q, units="gkg"):
    """spHumToMixRat(q)
    Calculate mixing ratio from specific humidity
    Assumes the input specific humidity variable is in units
    of g/kg.
    """
    q = convert(q, units, "kgkg")

    rat = gEps*q/(gEps-q)
    return rat

def rHToMixRat(rh, tmp, prs, tmp_units="C"):
    """rHToMixRat(rh, tmp, prs):
    Calculate mixing ratio from relative humidity, temperature and pressure
    """
    es = satVapPr(convert(tmp, tmp_units, "C"))
    e = (rh/100.)*es
    rat = vapPrToMixRat(e, prs)
    return rat

def spHumToRH(q, tmp, prs):
    """spHumToRH(tmp, prs):
    Calculate relative humidity from specific humidity, temperature and
    pressure.
    """
    es = satVapPr(tmp)
    qs = gEps * es / prs
    rh = 100.*q/qs
    return rh

def coriolis(lat):
    """Calculate the Coriolis factor
    Calculate the Coriolis factor (f) for a given latitude (degrees).
    If a list is passed, return a list, else return a single value.
    """
    omega = 2*math.pi/86400.
    f = 2*omega*numpy.sin(radians(lat))
    return f

def convert(value, input, output):
    """
    Convert value from input units to output units.
    """
    startValue = value
    value = numpy.array(value)
    if input == output:
        # Do nothing:
        return value
    if input=='kmh':
        input = 'kph'
    # Speeds:
    mps = {"kph":3.6, "kts":1.944, "mph":2.2369}
    mph = {"kph":1.60934, "kts":0.86898, "mps":0.44704}
    kph = {"kts":0.539957, "mps":0.2777778,"mph":0.621371}
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
    Pascals = {"kPa":0.001, "hPa":0.01, "inHg":0.0002953, "mmHg":0.007500616, 
          "Pa":1.0}

    # Lengths:
    km = {"m":0.001,"mi":0.621371192, "deg":0.00899886, "nm":0.539957,
          "rad":0.0001570783}
    deg = {"km":111.1251, "m":111125.1, "mi":69.0499358, "nm":60.0,
           "rad":math.pi/180.}
    mi = {"km":1.60934, "m":1609.34, "deg":0.014482}
    nm = {"km":1.852, "m":1852, "deg":0.01666, "rad":math.pi/10800.}
    rad = {"nm":10800./math.pi, "km":6366.248653, "deg":180./math.pi}

    # Mixing ratio:
    gkg = {"kgkg":0.001}
    kgkg = {"gkg":1000}

    convert={"mps":mps,
             "mph":mph,
             "kph":kph,
             "kts":kts,
             "kPa":kPa,
             "hPa":hPa,
             "Pa":Pa,
             "Pascals":Pascals,
             "inHg":inHg,
             "mmHg":mmHg,
             "C":C,
             "F":F,
             "K":K,
             "km":km,
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

    if input in convert_pre:
        if output in convert_pre[input]:
            value += convert_pre[input][output]

    if input in convert:
        if output in convert[input]:
            value = value*convert[input][output]

    if input in convert_post:
        if output in convert_post[input]:
            value += convert_post[input][output]

    return value

def vapour(temp):
    """
    Determine equivalent potential temperature (theta-e) given temperature
    """
    vapour = 6.112*numpy.exp(17.67*temp/(243.5+temp))
    return vapour

def genesisPotential(zeta, rh, vmax, shear):
    """
    Calculate genesis potential index
    """
    gpi = power(abs((10**5)*zeta),
                1.5)*((rh/50.)**3)*((vmax/70.)**3)/((1.+0.1*shear)**2)
    return gpi
