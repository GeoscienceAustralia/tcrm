subroutine beardist(cLon, cLat, lonArray, latArray, bearing, dist, nlon, nlat)
    !$ use omp_lib

!    :param double cLon: longitude of storm centre (degrees)
!    :param double clat: latitude of storm centre (degrees)
!    :param lonArray: 1D double precision array of longitudes (degrees)
!    :param latArray: 1D double precision array of latitudes (degrees)
!    :param bearing: 2D double precision output array of bearings with shape (nlat, nlon) (radians)
!    :param dist: 2D double precision output array of distances with shape (nlat, nlon) (km)
!    :param int nlon: length of lonArray
!    :param int nlat: length of latArray

    integer, intent(in) :: nlon, nlat
    doubleprecision, intent(in) :: lonArray(nlon), latArray(nlat)
    doubleprecision, intent(inout), dimension(nlat, nlon) :: bearing, dist

    doubleprecision :: toRads, cLon, cLat, dlon, lon(nlon), lat(nlat), radius
    doubleprecision :: dLon_sin(nlon), dLon_cos(nlon), lat_sin(nlat), lat_cos(nlat)
    doubleprecision :: dLat_sin(nlat), dhalfLon_sin(nlon), a, c, pi
    doubleprecision :: cLon_, cLat_, cLat_cos, cLat_sin, alpha, beta

    pi = 4.d0*datan(1.d0)
    toRads = 0.017453292519943295
    radius = 6367.0

    cLon_ = cLon;
    cLat_ = cLat;

    cLon_ = cLon_ * toRads
    cLat_ = cLat_ * toRads

    cLat_cos = cos(cLat_)
    cLat_sin = sin(cLat_)

    do i = 1, nlon
        lon(i) = lonArray(i) * toRads
        dLon = lon(i) - cLon_
        dLon_sin(i) = sin(dLon)
        dLon_cos(i) = cos(dLon)
        dhalfLon_sin(i) = sin(0.5 * dLon)
    end do

    do i = 1, nlat
        lat(i) = latArray(i) * toRads
        lat_sin(i) = sin(lat(i))
        lat_cos(i) = cos(lat(i))
        dLat_sin(i) = sin(0.5 * (lat(i) - cLat_))
    end do

    !$OMP PARALLEL DO shared(bearing, dist)
    do j = 1, nlat
        do i = 1, nlon

            alpha = dLon_sin(i) * lat_cos(j);
            beta = (cLat_cos * lat_sin(j)) - (cLat_sin * lat_cos(j) * dLon_cos(i));
            bearing(j, i) = 0.5 * pi - atan2(alpha, beta);

            a = dLat_sin(j) * dLat_sin(j) + cLat_cos * lat_cos(j) * dhalfLon_sin(i) * dhalfLon_sin(i);
            c = 2.0 * atan2(sqrt(abs(a)), sqrt(1.0 - a));
            dist(j, i) = max(radius*c, 1e-30);
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine beardist