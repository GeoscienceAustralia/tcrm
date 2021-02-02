"""
:mod:`lmomentFit` -- fit GEV functions
======================================

.. module:: lmomentFit
    :synopsis: Functions translated from the LMOMENTS Fortran
               package for fitting a GEV function.

.. moduleauthor:: Nicholas Summons <nicholas.summons@ga.gov.au>

Two functions {pelgev, samlmu} from the LMOMENTS Fortran package
ported to Python to fit a Generalised Extreme Value Distribution
function. Original code developed by: J. R. M. HOSKING, IBM RESEARCH
DIVISION, T. J. WATSON RESEARCH CENTER, YORKTOWN HEIGHTS, NEW YORK
10598, U.S.A.

.. note::
    Permission to use, copy, modify and distribute this software for
    any purpose and without fee is hereby granted, provided that this
    copyright and permission notice appear on all copies of the
    software. The name of the IBM Corporation may not be used in any
    advertising or publicity pertaining to the use of the
    software. IBM makes no warranty or representations about the
    suitability of the software for any purpose. It is provided "AS
    IS" without any express or implied warranty, including the implied
    warranties of merchantability, fitness for a particular purpose
    and non-infringement. IBM shall not be liable for any direct,
    indirect, special or consequential damages resulting from the loss
    of use, data or projects, whether in an action of contract or
    tort, arising out of or in connection with the use or performance
    of this software.

"""

import numpy as np
from scipy import special
import logging

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

__version__ = "$Id: lmomentFit.py 686 2012-03-29 04:24:59Z carthur $"

def pelgev(xmom):
    """
    Parameter estimation via L-moments for the Generalised Extreme
    Value Distribution. For -0.8 <= TAU3 < 1., K is approximated by
    rational functions as in Donaldson (1996,
    Commun. Statist. Simul. Comput.). If TAU3 is outside this range,
    Newton-Raphson iteration is used.

    :param xmom: Array of length 3, containing the L-moments Lambda-1,
                 Lambda-2 and TAU3.
    :type  xmom: List or :class:`numpy.ndarray`

    :returns: Location, scale and shape parameters of the GEV
              distribution.
    :rtype: :class:`numpy.ndarray`

    """
    #***********************************************************************
    #*                                                                     *
    #*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
    #*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
    #*                                                                     *
    #*  J. R. M. HOSKING                                                   *
    #*  IBM RESEARCH DIVISION                                              *
    #*  T. J. WATSON RESEARCH CENTER                                       *
    #*  YORKTOWN HEIGHTS                                                   *
    #*  NEW YORK 10598, U.S.A.                                             *
    #*                                                                     *
    #*  VERSION 3     AUGUST 1996                                          *
    #*                                                                     *
    #***********************************************************************
    #
    #  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED EXTREME-VALUE
    #  DISTRIBUTION
    #
    #  PARAMETERS OF ROUTINE:
    #  xmom   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
    #                  LAMBDA-2, TAU-3.
    #  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
    #                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
    #
    #  OTHER ROUTINES USED: DLGAMA
    #
    # METHOD: FOR  -0.8 LE TAU3 LT 1,  K IS APPROXIMATED BY RATIONAL
    # FUNCTIONS AS IN DONALDSON (1996, COMMUN. STATIST. SIMUL. COMPUT.).
    # IF TAU3 IS OUTSIDE THIS RANGE, NEWTON-RAPHSON ITERATION IS USED.
    #xmom = np.array([1.235, 0.11367, 0.10557])
    para = np.zeros(3)

    # small IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
    # epsilon, maxit CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
    small = 1E-5
    epsilon = 1E-6
    maxit = 20

    # EU IS EULER'S CONSTANT
    # DL2 IS LOG(2), DL3 IS LOG(3)

    EU = 0.57721566
    DL2 = 0.69314718
    DL3 = 1.0986123

    # COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS FOR K
    A0 = 0.28377530
    A1 = -1.21096399
    A2 = -2.50728214
    A3 = -1.13455566
    A4 = -0.07138022
    B1 = 2.06189696
    B2 = 1.31912239
    B3 = 0.25077104
    C1 = 1.59921491
    C2 = -0.48832213
    C3 = 0.01573152
    D1 = -0.64363929
    D2 = 0.08985247

    t3 = xmom[2]
    if xmom[1] <= 0.0:
        log.debug(' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID')
        return para
    if np.abs(t3) >= 1.0:
        log.debug(' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID')
        return para
    if t3 > 0.0:
        # RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN 0 AND 1
        Z = 1.0 - t3
        G = (-1.0 + Z * (C1 + Z * (C2 + Z * C3)))/(1.0 + Z * (D1 + Z * D2))
        if np.abs(G) < small:
            # ESTIMATED K EFFECTIVELY ZERO
            para[2] = 0.0
            para[1] = xmom[1] / DL2
            para[0] = xmom[0] - EU * para[1]
            return para
        para[2] = G
        gam = special.gamma(1.0 + G)
        para[1] = xmom[1] * G / (gam * (1.0 - 2.0**(-G)))
        para[0] = xmom[0] - para[1] * (1.0 - gam)/G
        return para

    # RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN -0.8 AND 0
    G = (A0 + t3 * (A1 + t3 * (A2 + t3 * (A3 + t3 * A4))))/\
        (1.0 + t3 * (B1 + t3 * (B2 + t3 * B3)))
    if t3 < -0.8:
        # NEWTON-RAPHSON ITERATION FOR TAU3 LESS THAN -0.8
        if t3 <= -0.97:
            G = 1.0 - np.log(1.0 + t3)/DL2
        T0 = (t3 + 3.0) * 0.5

        convg = False
        for IT in range(1, maxit + 1):
            x2 = 2.0**(-G)
            x3 = 3.0**(-G)
            xx2 = 1.0 - x2
            xx3 = 1.0 - x3
            T = xx3 / xx2
            DERIV = (xx2 * x3 * DL3 - xx3 * x2 * DL2)/(xx2 * xx2)
            GOLD = G
            G = G - (T - T0)/DERIV
            if np.abs(G - GOLD) <= epsilon * G:
                convg = True
                break

        if convg == False:
            log.debug((' ** WARNING ** ROUTINE PELGEV : ITERATION HAS '
                       'NOT CONVERGED. RESULTS MAY BE UNRELIABLE.'))

    # ESTIMATE ALPHA,XI
    para[2] = G
    gam = special.gamma(1.0 + G)
    para[1] = xmom[1] * G / (gam * (1.0 - 2.0**(-G)))
    para[0] = xmom[0] - para[1] * (1.0 - gam)/G
    return para

def pelgpa(xmom):
    """
    Parameter estimation via L-moments for the Generalised Pareto
    Distribution.
    """
    para = np.zeros(3)
    t3 = xmom[2]
    if xmom[1] <= 0.0:
        log.debug(' *** ERROR *** ROUTINE PELGPA : L-MOMENTS INVALID')
        return para
    if np.abs(t3) >= 1.0:
        log.debug(' *** ERROR *** ROUTINE PELGPA : L-MOMENTS INVALID')
        return para

    gg = (1.0 - 3.0 * t3) / (1.0 + t3)
    para[2] = gg
    para[1] = (1.0 + gg) * (2.0 + gg) * xmom[1]
    para[0] = xmom[0] - para[1] / (1.0 + gg)

    return para

def samlmu(data, nmom):
    """
    Sample L-moments for a data array.

    :param data: Array of length N, containing the data in ascending order.
    :type  data: :class:`numpy.ndarray`
    :param nmom: Number of L-moments to be found (maximum 100).

    :returns: The sample L-moments.
    :rtype: :class:`numpy.ndarray`
    """
    #***********************************************************************
    #*                                                                     *
    #*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
    #*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
    #*                                                                     *
    #*  J. R. M. HOSKING                                                   *
    #*  IBM RESEARCH DIVISION                                              *
    #*  T. J. WATSON RESEARCH CENTER                                       *
    #*  YORKTOWN HEIGHTS                                                   *
    #*  NEW YORK 10598, U.S.A.                                             *
    #*                                                                     *
    #*  VERSION 3     AUGUST 1996                                          *
    #*                                                                     *
    #***********************************************************************
    #
    #  SAMPLE L-MOMENTS OF A DATA ARRAY
    #
    #  PARAMETERS OF ROUTINE:
    #  x      * INPUT* ARRAY OF LENGTH N. CONTAINS THE DATA, IN ASCENDING
    #                  ORDER.
    #  N      * INPUT* NUMBER OF DATA VALUES
    #  xmom   *OUTPUT* ARRAY OF LENGTH nmom. CONTAINS THE SAMPLE L-MOMENTS,
    #                  STORED AS DESCRIBED BELOW.
    #  nmom   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 100.
    #

    if len(data) <= 2:
        raise ValueError("Not enough values for L-moment calculation")
    # If nmom == 3, use optimised code
    if nmom == 3:
        return samlmu3(data)

    maxmom = 100
    data = np.array(data)
    n = np.size(data)
    nmom = int(nmom)
    xmom = np.zeros(nmom)
    coef = np.zeros([2, nmom])

    if nmom > maxmom:
        log.debug(' ** WARNING ** ROUTINE SAMLMU : PARAMETER nmom INVALID')
        return
    dn = n
    xmom[:] = 0.0
    if nmom <= 2.0:
        # AT MOST TWO L-MOMENTS
        sum1 = 0.0
        sum2 = 0.0
        temp = -dn + 1.0
        for i in range(1, n + 1):
            sum1 = sum1 + data[i - 1]
            sum2 = sum2 + data[i - 1] * temp
            temp = temp + 2.0
        xmom[1-1] = sum1 / dn
        if nmom == 1:
            return xmom
        xmom[1] = sum2 / (dn * (dn - 1.0))
        return xmom

    # UNBIASED ESTIMATES OF L-MOMENTS -- THE 'DO 30' LOOP
    # RECURSIVELY CALCULATES DISCRETE LEGENDRE POLYNOMIALS, VIA
    # EQ.(9) OF NEUMAN AND SCHONBACH (1974, INT.J.NUM.METH.ENG.)
    for j in range(3, nmom + 1):
        temp = 1.0/((j - 1) * (n - j + 1))
        coef[0, j-1] = (j + j - 3) * temp
        coef[1, j-1] = ((j - 2) * (n + j - 2)) * temp

    temp = -dn - 1.0
    const = 1.0 / (dn - 1.0)
    nhalf = n//2
    for i in range(1, nhalf + 1):
        temp = temp + 2.0
        xi = data[i - 1]
        xii = data[n - i]
        termp = xi + xii
        termn = xi - xii
        xmom[0] = xmom[0] + termp
        S1 = 1.0
        S = temp * const
        xmom[1] = xmom[1] + S * termn
        for j in range(3, nmom + 1, 2):
            S2 = S1
            S1 = S
            S = coef[0, j - 1] * temp * S1 - coef[1, j - 1] * S2
            xmom[j - 1] = xmom[j - 1] + S * termp
            if j == nmom:
                break
            jj = j + 1
            S2 = S1
            S1 = S
            S = coef[0, jj - 1] * temp * S1 - coef[1, jj - 1] * S2
            xmom[jj - 1] = xmom[jj - 1] + S * termn

    if not (n == nhalf + nhalf):
        term = data[nhalf]
        s = 1.0
        xmom[0] = xmom[0] + term
        for j in range(3, nmom + 1, 2):
            s = -coef[1, j - 1] * s
            xmom[j - 1] = xmom[j - 1] + s * term

    # L-MOMENT RATIOS
    xmom[0] = xmom[0]/dn
    if xmom[1] == 0.0:
        log.debug(' *** ERROR *** ROUTINE SAMLMU : ALL DATA VALUES EQUAL')
        xmom[:] = 0.0
        return
    for j in range(3, nmom + 1):
        xmom[j - 1] = xmom[j - 1] / xmom[1]
    xmom[1] = xmom[1] / dn
    return np.array(xmom)


def samlmu3(data):
    """
    Functional equivalent to lmoments.samlmu(data, 3). Vectorised for
    speed but still about half speed of original fortran package.

    :param data: Array of length N, containing the data in ascending order.
    :returns: First 3 L-moments.
    :rtype: :class:`numpy.ndarray`

    """
    if len(data) <= 2:
        raise ValueError("Not enough values for L-moment calculation")
    data = np.array(data)
    n = np.size(data)
    xmom = np.zeros(3)

    tempo = 1.0 / (2 * n - 4)
    coef02 = 3 * tempo
    coef12 = (n + 1) * tempo

    temp = list(range(-n + 1, 0, 2))
    const = 1.0 / (n - 1.0)
    nhalf = n // 2
    xi = data[0:nhalf]
    xii = data[-1:-nhalf - 1:-1]
    termp = xi + xii
    termn = xi - xii
    xmom[0] = sum(termp)
    xmom[1] = sum(temp * termn * const)
    s = [k * k * coef02 * const - coef12 for k in temp]
    xmom[2] = sum(s * termp)

    if n != nhalf + nhalf:
        xmom[0] = xmom[0] + data[nhalf]
        xmom[2] = xmom[2] - coef12 * data[nhalf]
    xmom = [xmom[0] / n, xmom[1] / n, xmom[2] / xmom[1]]
    return np.array(xmom)

