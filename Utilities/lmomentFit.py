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

import numpy
from scipy import special
__version__ = "$Id: lmomentFit.py 686 2012-03-29 04:24:59Z carthur $"

def pelgev(XMOM):
    """
    Parameter estimation via L-moments for the Generalised Extreme
    Value Distribution. For -0.8 <= TAU3 < 1., K is approximated by
    rational functions as in Donaldson (1996,
    Commun. Statist. Simul. Comput.). If TAU3 is outside this range,
    Newton-Raphson iteration is used.

    :param XMOM: Array of length 3, containing the L-moments Lambda-1,
                 Lambda-2 and TAU3.
    :type  XMOM: List or :class:`numpy.ndarray`

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
    #  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
    #                  LAMBDA-2, TAU-3.
    #  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
    #                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
    #
    #  OTHER ROUTINES USED: DLGAMA
    #
    # METHOD: FOR  -0.8 LE TAU3 LT 1,  K IS APPROXIMATED BY RATIONAL
    # FUNCTIONS AS IN DONALDSON (1996, COMMUN. STATIST. SIMUL. COMPUT.).
    # IF TAU3 IS OUTSIDE THIS RANGE, NEWTON-RAPHSON ITERATION IS USED.
    #XMOM = numpy.array([1.235, 0.11367, 0.10557])
    PARA = numpy.zeros(3)
    P8 = 0.8
    P97 = 0.97

    # SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
    # EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
    SMALL = 1E-5
    EPS = 1E-6
    MAXIT = 20

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

    T3 = XMOM[2]
    if XMOM[1] <= 0.0:
        print ' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID'
        return PARA
    if numpy.abs(T3) >= 1.0:
        print ' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID'
        return PARA
    if T3 > 0.0:
        # RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN 0 AND 1
        Z = 1.0 - T3
        G = (-1.0+Z*(C1+Z*(C2+Z*C3)))/(1.0+Z*(D1+Z*D2))
        if numpy.abs(G) < SMALL:
            # ESTIMATED K EFFECTIVELY ZERO
            PARA[2] = 0.0
            PARA[1] = XMOM[1]/DL2
            PARA[0] = XMOM[0] - EU*PARA[1]
            return PARA
        PARA[2] = G
        GAM = special.gamma(1.0+G)
        PARA[1] = XMOM[1]*G/(GAM*(1.0-2.0**(-G)))
        PARA[0] = XMOM[0]-PARA[1]*(1.0-GAM)/G
        return PARA

    # RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN -0.8 AND 0
    G = (A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(1.0+T3*(B1+T3*(B2+T3*B3)))
    if T3 < -P8:
        # NEWTON-RAPHSON ITERATION FOR TAU3 LESS THAN -0.8
        if T3 <= -P97:
            G = 1.0 - numpy.log(1.0 + T3)/DL2
        T0 = (T3 + 3.0)*0.5

        convg = False
        for IT in xrange(1, MAXIT+1):
            X2 = 2.0**(-G)
            X3 = 3.0**(-G)
            XX2 = 1.0 - X2
            XX3 = 1.0 - X3
            T = XX3/XX2
            DERIV = (XX2*X3*DL3 - XX3*X2*DL2)/(XX2*XX2)
            GOLD = G
            G = G - (T - T0)/DERIV
            if numpy.abs(G-GOLD) <= EPS*G:
                convg = True
                break

        if convg == False:
            print ' ** WARNING ** ROUTINE PELGEV : ITERATION HAS NOT CONVERGED. RESULTS MAY BE UNRELIABLE.'

    # ESTIMATE ALPHA,XI
    PARA[2] = G
    GAM = special.gamma(1.0+G)
    PARA[1] = XMOM[1]*G/(GAM*(1.0-2.0**(-G)))
    PARA[0] = XMOM[0]-PARA[1]*(1.0-GAM)/G
    return PARA

def pelgpa(XMOM):
    
    PARA = numpy.zeros(3)
    T3 = XMOM[2]
    if XMOM[1] <= 0.0:
        print ' *** ERROR *** ROUTINE PELGPA : L-MOMENTS INVALID'
        return PARA
    if numpy.abs(T3) >= 1.0:
        print ' *** ERROR *** ROUTINE PELGPA : L-MOMENTS INVALID'
        return PARA

    G = (1.0 - 3.0 * T3) / (1.0 + T3)
    PARA[2] = G
    PARA[1] = (1.0 + G) * (2.0 + G) * XMOM[1]
    PARA[0] = XMOM[0] - PARA[1] / (1.0 + G)

    return PARA

def samlmu(X, NMOM):
    """
    Sample L-moments for a data array.

    :param X: Array of length N, containing the data in ascending order.
    :type  X: :class:`numpy.ndarray`
    :param NMOM: Number of L-moments to be found (maximum 100).

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
    #  X      * INPUT* ARRAY OF LENGTH N. CONTAINS THE DATA, IN ASCENDING
    #                  ORDER.
    #  N      * INPUT* NUMBER OF DATA VALUES
    #  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. CONTAINS THE SAMPLE L-MOMENTS,
    #                  STORED AS DESCRIBED BELOW.
    #  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 100.
    #

    # If NMOM == 3, use optimised code
    if NMOM == 3:
        return samlmu3(X)

    MAXMOM = 100
    X = numpy.array(X)
    N = numpy.size(X)
    NMOM = int(NMOM)
    XMOM = numpy.zeros(NMOM)
    COEF = numpy.zeros([2,NMOM])

    if NMOM > MAXMOM:
        print ' ** WARNING ** ROUTINE SAMLMU : PARAMETER NMOM INVALID'
        return
    DN = N
    XMOM[:] = 0.0
    if NMOM <= 2.0:
        # AT MOST TWO L-MOMENTS
        SUM1 = 0.0
        SUM2 = 0.0
        TEMP = -DN + 1.0
        for I in xrange(1,N+1):
            SUM1 = SUM1 + X[I-1]
            SUM2 = SUM2 + X[I-1]*TEMP
            TEMP = TEMP + 2.0
        XMOM[1-1] = SUM1/DN
        if NMOM == 1:
            return XMOM
        XMOM[1] = SUM2/(DN*(DN-1.0))
        return XMOM

    # UNBIASED ESTIMATES OF L-MOMENTS -- THE 'DO 30' LOOP
    # RECURSIVELY CALCULATES DISCRETE LEGENDRE POLYNOMIALS, VIA
    # EQ.(9) OF NEUMAN AND SCHONBACH (1974, INT.J.NUM.METH.ENG.)
    for J in xrange(3, NMOM + 1):
        TEMP = 1.0/((J - 1)*(N - J + 1))
        COEF[0,J-1] = (J + J - 3)*TEMP
        COEF[1,J-1] = ((J - 2)*(N + J - 2))*TEMP

    TEMP = -DN - 1.0
    CONST = 1.0/(DN - 1.0)
    NHALF = N/2
    for I in xrange(1,NHALF + 1):
        TEMP = TEMP + 2.0
        XI = X[I-1]
        XII = X[N - I]
        TERMP = XI + XII
        TERMN = XI - XII
        XMOM[0] = XMOM[0] + TERMP
        S1 = 1.0
        S = TEMP*CONST
        XMOM[1] = XMOM[1] + S*TERMN
        for J in xrange(3,NMOM + 1,2):
            S2 = S1
            S1 = S
            S = COEF[0,J - 1]*TEMP*S1 - COEF[1,J - 1]*S2
            XMOM[J - 1] = XMOM[J - 1] + S*TERMP
            if J == NMOM:
                break
            JJ = J + 1
            S2 = S1
            S1 = S
            S = COEF[0,JJ - 1]*TEMP*S1 - COEF[1,JJ - 1]*S2
            XMOM[JJ - 1] = XMOM[JJ - 1] + S*TERMN

    if not (N == NHALF+NHALF):
        TERM = X[NHALF]
        S = 1.0
        XMOM[0] = XMOM[0] + TERM
        for J in xrange(3,NMOM+1,2):
            S = -COEF[1, J - 1]*S
            XMOM[J - 1] = XMOM[J - 1] + S*TERM

    # L-MOMENT RATIOS
    XMOM[0] = XMOM[0]/DN
    if XMOM[1] == 0.0:
        print ' *** ERROR *** ROUTINE SAMLMU : ALL DATA VALUES EQUAL'
        XMOM[:] = 0.0
        return
    for J in xrange(3, NMOM + 1):
        XMOM[J - 1] = XMOM[J - 1]/XMOM[1]
    XMOM[1] = XMOM[1]/DN
    return numpy.array(XMOM)


def samlmu3(X):
    """
    Functional equivalent to lmoments.samlmu(X, 3). Vectorised for
    speed but still about half speed of original fortran package.

    :param X: Array of length N, containing the data in ascending order.
    :returns: First 3 L-moments.
    :rtype: :class:`numpy.ndarray`
    
    """
    X = numpy.array(X)
    N = numpy.size(X)
    XMOM = numpy.zeros(3)

    TEMP0 = 1.0/(2*N - 4)
    COEF02 = 3*TEMP0
    COEF12 = (N + 1)*TEMP0

    TEMP = range(-N+1,0,2)
    CONST = 1.0/(N - 1.0)
    NHALF = N/2
    XI = X[0:NHALF]
    XII = X[-1:-NHALF-1:-1]
    TERMP = XI + XII
    TERMN = XI - XII
    XMOM[0] = sum(TERMP)
    XMOM[1] = sum(TEMP*TERMN*CONST)
    S = [k*k*COEF02*CONST - COEF12 for k in TEMP]
    XMOM[2] = sum(S*TERMP)

    if N != NHALF+NHALF:
        XMOM[0] = XMOM[0] + X[NHALF]
        XMOM[2] = XMOM[2] - COEF12*X[NHALF]
    XMOM = [XMOM[0]/N, XMOM[1]/N, XMOM[2]/XMOM[1]]
    return numpy.array(XMOM)
