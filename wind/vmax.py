"""
:mod:`vmax` -- calculate maximum wind speed from pressure deficit
=================================================================

.. module:: vmax
    :synopsis: Calculate maximum wind speed in a cyclone from the
               pressure deficit.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>


CreationDate: 2006-10-30
Description: Calculate Vmax from pressure difference and vice versa.
Offers three relations for calculating Vmax based on the pressure
difference -
Atkinson & Holliday (1977), Holland (1980) and Willoughby & Rahn (2004).
The default is the Holland relation.

The normalisations to 10m AGL, 10-minute-mean winds were taken from
Harper (2002).

SeeAlso: (related programs)
Constraints: Input pressures should be in Pa (not hPa)
Version: $Rev: 810 $

Version: 220
ModifiedBy: Craig Arthur
ModifiedDate: 2006-12-01
Modification: Added pDiff function

Version: 267
ModifiedBy: Craig Arthur
ModifiedDate: 2007-02-19
Modification: Normalised all vmax relations to provide 10m AGL,
10-minute-mean wind speeds

References
 * Atkinson and Holliday 1977: Tropical Cyclone Minimum Sea Level Pressure /
   Maximum Sustained Wind Relationship for the Wetern North Pacific.
   Mon. Wea. Rev., 105, 421-427
 * Holland, G. 1980: An Analytic Model of the Wind and Pressure Profiles in
   Hurricanes.  Mon. Wea. Rev, 108, 1212-1218
 * Willoughby, H.E. and M.E. Rahn, 2004: Parametric Representation of the
   Primary Hurricane Vortex. Part I: Observations and Evaluation of the Holland
   (1980) Model. Mon. Wea. Rev., 132, 3033-3048
 * Harper, B.A. 2002: Tropical Cyclone Parameter Estimation in the Australian
   Region: Wind-Pressure Relationships and Related Issues for Engineering
   Planning and Design.

$Id: vmax.py 810 2012-02-21 07:52:50Z nsummons $

"""

from scipy import sqrt, exp, power
import Utilities.metutils as metutils


def vmax(pCentre, pEnv, type="holland", beta=1.3, rho=1.15):
    """
    Calculate the maximum wind speed from the pressure difference.

    :param float pc: central pressure (Pa)
    :param float pe: environmental pressure (Pa)
    :param str type: which Vmax relation to use (Willoughby & Rahn,
                     Holland or Atkinson & Holliday)
    :param float beta: Holland's (1980) beta parameter. Only used for the
                       Holland estimation (type=holland)
    :param float rho: air density (default=1.15 kg/m^3)
    :return: maximum wind speed. For types 1 & 2, this is a
             gradient level wind. The relation used in type 3
             (Atkinson & Holliday) was determined using surface
             wind observations so should be used with caution at
             the gradient level.
    :raises ValueError: if environmental pressure is lower than
                        central pressure

    Note: The pressure should ideally be passed in units of Pa, but the
    function will accept hPa and automatically convert to Pa.

    """
    # Convert from hPa to Pa if necessary:
    if pCentre < 10000:
        pCentre = metutils.convert(pCentre, "hPa", "Pa")

    if pEnv < 10000:
        pEnv = metutils.convert(pEnv, "hPa", "Pa")

    if pEnv < pCentre:
        raise ValueError(("Error in vmax - Environmental pressure is less "
                          "than central pressure. Check values and/or order "
                          "of input arguments"))

    dP = pEnv - pCentre

    if type == "willoughby":
        # Default: Most advanced estimation technique:
        # Willoughby & Rahn (2004), Parametric Representation of the
        # Primary Hurricane Vortex. Part I: Observations and
        # Evaluation of the Holland (1980) Model.
        # Mon. Wea. Rev., 132, 3033-3048
        vMax = 0.6252*sqrt(dP)
    elif type == "holland":
        # Holland (1980), An Analytic Model of the Wind and Pressure
        # Profiles in Hurricanes. Mon. Wea. Rev, 108, 1212-1218
        # Density of air is assumed to be 1.15 kg/m^3.
        # beta is assumed to be 1.3. Other values can be specified.
        # Gradient level wind (assumed maximum).
        vMax = sqrt(beta*dP/(exp(1)*rho))
    elif type == "atkinson":
        # Atkinson and Holliday (1977), Tropical Cyclone Minimum Sea
        # Level Pressure / Maximum Sustained Wind Relationship for
        # the Western North Pacific. Mon. Wea. Rev., 105, 421-427
        # Maximum 10m, 1-minute wind speed. Uses pEnv as 1010 hPa
        vMax = 3.04*power(1010 - metutils.convert(pCentre, "Pa", "hPa"), 0.644)
    else:
        raise NotImplementedError("Vmax type " + type + " not implemented")
    return vMax


def pDiff(vMax, pEnv, vMaxType="holland", beta=1.3, rho=1.15):
    """
    Inverse functions to calculate central pressure from vMax
    Assumes vMax is given in metres/second.
    Returns pCentre in Pa.

    """
    if pEnv < 10000:
        pEnv = metutils.convert(pEnv, "hPa", "Pa")

    if vMaxType == "willoughby":
        dP = (vMax/0.6252)**2
    elif vMaxType == "holland":
        dP = rho*exp(1)*(vMax**2)/beta
    elif vMaxType == "atkinson":
        dP = (vMax/3.04)**(1/0.644)
        dP = metutils.convert(dP, "hPa", "Pa")
    else:
        raise NotImplementedError("Vmax type " + vMaxType + " not implemented")

    pCentre = pEnv - dP
    return pCentre
