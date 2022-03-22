subroutine fkerpert(R, lam, V, Z, f, rMax, vFm, thetaFm, Vm, Ux, Uy, n)
   !$ use omp_lib

   integer, intent(in) :: n
   doubleprecision, intent(in) :: f, rMax, vFm, Vm
   doubleprecision, dimension(n), intent(in) :: R, lam, V, Z
   doubleprecision, dimension(n), intent(inout) :: Ux, Uy

   doubleprecision :: Umod, Vt, al, be, gam , K, Cd, u0s, v0s, chi, eta, psi, albe, b
   doubleprecision :: ups, ums, vps, vms, usf, vsf, phi, us, vs
   complex(8) :: A0, j, Am, Ap
   logical :: ind

   b = 1.0
   K = 50.0
   Cd = 0.002
   j = cmplx(0.0, 1.0)

   !$OMP PARALLEL DO shared(Ux, Uy)
   do i = 1, n

      if ((vFm > 0) .and. (Vm/vFm < 5.)) then
         Umod = vFm * abs(1.25*(1. - (vFm/Vm)))
      else
         Umod = vFm
      end if

      if (R(i) > 2 * rMax) then
         Vt  = Umod * exp(-((R(i) / (2.*rMax)) - 1.) ** 2.)
      else
         Vt = Umod
      end if

      al = ((2. * V(i) / R(i)) + f) / (2. * K)
      be = (f + Z(i)) / (2. * K)
      gam = (V(i) / (2. * K * R(i)))

      albe = sqrt(al / be)

      ind = abs(gam) > sqrt(al * be)
      chi = abs((Cd / K) * V(i) / sqrt(sqrt(al * be)))
      eta = abs((Cd / K) * V(i) / sqrt(sqrt(al * be) + abs(gam)))
      psi = abs((Cd / K) * V(i) / sqrt(abs(sqrt(al * be) - abs(gam))))

      A0 = -(chi * (1.0 + j * (1.0 + chi)) * V(i)) / (2.0 * chi**2 + 3.0 * chi + 2.0)

      u0s = realpart(A0) * albe * sign(b, f)
      v0s = imagpart(A0)

      if (ind) then
         Am = -(psi * (1 + 2 * albe + (1 + j) * (1 + albe) * eta) * Vt)
         Am = Am / (albe * ((2 - 2 * j + 3 * (eta + psi) + (2 + 2 * j) * eta * psi)))
         Ap = -(eta * (1 - 2 * albe + (1 - j) * (1 - albe)*psi) * Vt)
         Ap = Ap / (albe * (2 + 2 * j + 3 * (eta + psi) + (2 - 2 * j) * eta * psi))
      else
         Am = -(psi * (1 + 2 * albe + (1 + j) * (1 + albe) * eta) * Vt)
         Am = Am / (albe * ((2 + 2 * j) * (1 + eta * psi) + 3 * psi + 3 * j * eta))
         Ap = -(eta * (1 - 2 * albe + (1 + j) * (1 - albe) * psi) * Vt)
         Ap = Ap / (albe * ((2 + 2 * j) * (1 + eta * psi) + 3 * eta + 3 * j * psi))
      end if

      ! First asymmetric surface component
      ums = realpart(Am * exp(-j * (lam(i) - thetaFm) * sign(b, f))) * albe
      vms = imagpart(Am * exp(-j * (lam(i) - thetaFm) * sign(b, f))) * sign(b, f)

      ! Second asymmetric surface component
      ups = realpart(Ap * exp(j * (lam(i) - thetaFm) * sign(b, f))) * albe
      vps = realpart(Ap * exp(j * (lam(i) - thetaFm) * sign(b, f))) * sign(b, f)

      ! Total surface wind in (moving coordinate system)
      us = u0s + ups + ums
      vs = v0s + vps + vms + V(i)

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