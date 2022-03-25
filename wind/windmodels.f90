subroutine fkerpert(R, lam, f, rMax, Vm, thetaFm, vFm, d2Vm, dVm, dP, beta, rho, Ux, Uy, n)
   !$ use omp_lib

!    Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the
!    Tropical Cyclone Core. Part I: Linear Theory.  J. Atmos. Sci., 58,
!    2469-2484

!    This calculates the Kepert wind field with a Holland pressure profile.

!   :param R: Distance from the storm centre to the grid (m)
!     :type  R: 1D double precision array
!   :param lam: Direction (0=east, radians, positive anti-clockwise) from storm centre to the grid.
!     :type  lam: 1D double precision array
!   :param float rMax: Radius to maximum gradient winds (m).
!   :param float Vm: Maximum gradient wind speed (m/s).
!   :param float thetaFm: Bearing of storm (0=east, radians positive anti-clockwise)..
!   :param float vFm: Foward speed of the storm (m/s).
!   :param float d2Vm: Second derivitive of gradient windpseed w.r.t. radius at rMax.
!   :param float dVm: Derivitive of gradient windpseed w.r.t. radius at rMax.
!   :param float dP: Central pressure deficit in Pa
!   :param float beta: Holland beta parameter.
!   :param float rho: Density of air (kg/m^3).
!   :param Ux: Output east windspeed (m/s)
!     :type  Ux: 1D double precision array
!   :param Uy: Output north windspeed (m/s)
!     :type  Uy: 1D double precision array
!   :param n: length of arrays

   integer, intent(in) :: n
   doubleprecision, intent(in) :: f, rMax, Vm, thetaFm, vFm, d2Vm, dVm, dP, beta, rho
   doubleprecision, dimension(n), intent(in) :: R, lam
   doubleprecision, dimension(n), intent(inout) :: Ux, Uy

   doubleprecision :: Umod, Vt, al, be, gam , K, Cd, u0s, v0s, chi, eta, psi, albe, b
   doubleprecision :: ups, ums, vps, vms, usf, vsf, phi, us, vs
   complex(8) :: A0, j, Am, Ap
   logical :: ind
   doubleprecision :: aa, bb, cc, delta, edelta, V, Z
   logical :: icore

   b = 1.0
   K = 50.0
   Cd = 0.002
   j = cmplx(0.0, 1.0, 8)

   aa = (d2Vm / 2 - (dVm - Vm / rMax) / rMax) / rMax
   bb = (d2Vm - 6 * aa * rMax) / 2.
   cc = dVm - 3 * aa * rMax ** 2 - 2 * bb * rMax


   if ((vFm > 0) .and. (Vm/vFm < 5.)) then
         Umod = vFm * abs(1.25*(1. - (vFm/Vm)))
   else
      Umod = vFm
   end if

   !$OMP PARALLEL DO shared(Ux, Uy)
   do i = 1, n

      delta = (rMax / R(i)) ** beta
      edelta = exp(-delta)
      icore = (R(i) <= rMax)
      if (icore) then
         Z = R(i) * (R(i) * 4 * aa + 3 * bb) + 2 * cc
         V = (R(i) * (R(i) * (R(i) * aa + bb) + cc))
      else
         V = sqrt((dP * beta / rho) * delta * edelta + (R(i) * f / 2.) ** 2) - R(i) * abs(f) / 2.
         Z = abs(f) + &
            (beta**2 * dP * (delta**2) * edelta / &
             (2 * rho * R(i)) - beta**2 * dP * delta * edelta / &
             (2 * rho * R(i)) + R(i) * f**2 / 4) / &
            sqrt(beta * dP * delta * edelta / &
                    rho + (R(i) * f / 2)**2) + &
            (sqrt(beta * dP * delta * edelta / &
                     rho + (R(i) * f / 2)**2)) / R(i)
      end if
      V = sign(V, f)
      Z = sign(Z, f)

      if (R(i) > 2 * rMax) then
         Vt  = Umod * exp(-((R(i) / (2.*rMax)) - 1.) ** 2.)
      else
         Vt = Umod
      end if

      al = ((2. * V / R(i)) + f) / (2. * K)
      be = (f + Z) / (2. * K)
      gam = (V / (2. * K * R(i)))

      albe = sqrt(al / be)

      ind = abs(gam) > sqrt(al * be)
      chi = abs((Cd / K) * V / sqrt(sqrt(al * be)))
      eta = abs((Cd / K) * V / sqrt(sqrt(al * be) + abs(gam)))
      psi = abs((Cd / K) * V / sqrt(abs(sqrt(al * be) - abs(gam))))

      A0 = -(chi * (1.0 + j * (1.0 + chi)) * V) / (2.0 * chi**2 + 3.0 * chi + 2.0)

      u0s = realpart(A0) * albe * sign(b, f)
      v0s = imagpart(A0)

      if (ind) then
         Am = -(psi * (1 + 2 * albe + (1 + j) * (1 + albe) * eta) * Vt)
         Am = Am / (albe * ((2 - 2 * j + 3 * (eta + psi) + (2 + 2 * j) * eta * psi)))
         Ap = -(eta * (1 - 2 * albe + (1 - j) * (1 - albe) * psi) * Vt)
         Ap = Ap / (albe * (2 + 2 * j + 3 * (eta + psi) + (2 - 2 * j) * eta * psi))
      else
         Am = -(psi * (1 + 2 * albe + (1 + j) * (1 + albe) * eta) * Vt)
         Am = Am / (albe * ((2 + 2 * j) * (1 + eta * psi) + 3 * psi + 3 * j * eta))
         Ap = -(eta * (1.0 - 2.0 * albe + (1.0 + j) * (1.0 - albe) * psi) * Vt)
         Ap = Ap / (albe * ((2 + 2 * j) * (1 + eta * psi) + 3 * eta + 3 * j * psi))
      end if

      ! First asymmetric surface component
      ums = realpart(Am * exp(-j * (lam(i) - thetaFm) * sign(b, f))) * albe
      vms = imagpart(Am * exp(-j * (lam(i) - thetaFm) * sign(b, f))) * sign(b, f)

      ! Second asymmetric surface component
      ups = realpart(Ap * exp(j * (lam(i) - thetaFm) * sign(b, f))) * albe
      vps = imagpart(Ap * exp(j * (lam(i) - thetaFm) * sign(b, f))) * sign(b, f)

      ! Total surface wind in (moving coordinate system)
      us = u0s + ups + ums
      vs = v0s + vps + vms + V

      usf = us + Vt * cos(lam(i) - thetaFm)
      vsf = vs - Vt * sin(lam(i) - thetaFm)
      phi = atan2(usf, vsf)

      ! Surface winds, cartesian coordinates
      Ux(i) = sqrt(usf ** 2. + vsf ** 2.) * sin(phi - lam(i))
      Uy(i) = sqrt(usf ** 2. + vsf ** 2.) * cos(phi - lam(i))

   end do
   !$OMP END PARALLEL DO

end subroutine fkerpert


subroutine fhollandvel(V, R, d2Vm, dVm, rMax, vMax, beta, dP, rho, f, n)
   !$ use omp_lib
   doubleprecision, intent(in) :: d2Vm, dVm, rMax, beta, dP, rho, f, vMax
   integer, intent(in) :: n
   doubleprecision, intent(in), dimension(n) :: R
   doubleprecision, intent(inout), dimension(n) :: V
   doubleprecision :: aa, bb, cc, delta, edelta
   logical :: icore

   aa = ((d2Vm / 2. - (dVm - vMax / rMax) / rMax) / rMax)
   bb = (d2Vm - 6 * aa * rMax) / 2.
   cc = dVm - 3 * aa * rMax ** 2 - 2 * bb * rMax

   !$OMP PARALLEL DO shared(V)
   do i = 1, n

      delta = (rMax / R(i)) ** beta
      edelta = exp(-delta)

      icore = R(i) <= rMax
      if (icore) then
         V(i) = (R(i) * (R(i) * (R(i) * aa + bb) + cc))
      else
         V(i) = sqrt((dP * beta / rho) * delta * edelta + (R(i) * f / 2.) ** 2) - R(i) * abs(f) / 2.
      end if

      V(i) = sign(V(i), f)

   end do
   !$OMP END PARALLEL DO

end subroutine fhollandvel


subroutine fhollandvort(Z, R, d2Vm, dVm, rMax, vMax, beta, dP, rho, f, n)
   !$ use omp_lib
   doubleprecision, intent(in) :: d2Vm, dVm, rMax, beta, dP, rho, f, vMax
   integer, intent(in) :: n
   doubleprecision, intent(in), dimension(n) :: R
   doubleprecision, intent(inout), dimension(n) :: Z
   doubleprecision :: aa, bb, cc, delta, edelta
   logical :: icore

   aa = (d2Vm / 2 - (dVm - vMax / rMax) / rMax) / rMax
   bb = (d2Vm - 6 * aa * rMax) / 2.
   cc = dVm - 3 * aa * rMax ** 2 - 2 * bb * rMax

   !$OMP PARALLEL DO shared(Z)
   do i = 1, n
      delta = (rMax / R(i)) ** beta
      edelta = exp(-delta)
      icore = (R(i) <= rMax)
      if (icore) then
         Z(i) = R(i) * (R(i) * 4 * aa + 3 * bb) + 2 * cc
      else
         Z(i) = abs(f) + &
            (beta**2 * dP * (delta**2) * edelta / &
             (2 * rho * R(i)) - beta**2 * dP * delta * edelta / &
             (2 * rho * R(i)) + R(i) * f**2 / 4) / &
            sqrt(beta * dP * delta * edelta / &
                    rho + (R(i) * f / 2)**2) + &
            (sqrt(beta * dP * delta * edelta / &
                     rho + (R(i) * f / 2)**2)) / R(i)
      end if

      Z(i) = sign(Z(i), f)
   end do
   !$OMP END PARALLEL DO

end subroutine fhollandvort
