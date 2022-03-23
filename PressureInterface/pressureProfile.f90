subroutine fhollandpressure(P, R, rMax, pc, dP, beta, n)
   !$ use omp_lib
   integer, intent(in) :: n
   doubleprecision, intent(in), dimension(n) :: R
   doubleprecision, intent(inout), dimension(n) :: P
   doubleprecision, intent(in) :: rMax, beta, pc, dP

   !$OMP PARALLEL DO shared(P)
   do i = 1, n
      P(i) = pCentre + dP * exp(-(rMax / R(i)) ** beta)
   end do
   !$OMP END PARALLEL DO

end subroutine fhollandpressure