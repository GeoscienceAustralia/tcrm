subroutine fhollandpressure(P, R, rMax, pc, dP, beta, n)
   !$ use omp_lib

!   Calculates the pressure of the Holland profile
!
!   :param P: 1D double precision output pressure array
!   :param R: 1D double precision pressure array of distance from storm centre
!   :param double rMax: radius to maximum winds (m)
!   :param double pc: central pressure (Pa)
!   :param double dP: central pressure deficit (Pa)
!   :param double beta: shape parameter
!   :param int n: length of arrays

   integer, intent(in) :: n
   doubleprecision, intent(in), dimension(n) :: R
   doubleprecision, intent(inout), dimension(n) :: P
   doubleprecision, intent(in) :: rMax, beta, pc, dP

   !$OMP PARALLEL DO shared(P)
   do i = 1, n
      P(i) = pc + dP * exp(-(rMax / R(i)) ** beta)
   end do
   !$OMP END PARALLEL DO

end subroutine fhollandpressure