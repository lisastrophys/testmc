program part1
implicit none
integer, PARAMETER :: DP = KIND(1.0d0)
real(DP), parameter :: RSun   = 6.957e10_DP                 !Solar radius              [cm]
real(DP), parameter :: clight = 2.99792458e10_DP            !Speed of light            [cm s-1]
real(DP), parameter :: hP     = 6.626070040e-27_DP          !Planck constant           [erg s]
real(DP), parameter :: kB     = 1.38064852e-16_DP           !Boltzmann constant        [erg K-1]
real(DP), parameter :: AU     = 1.495978707e13_DP           !Astronomical unit         [cm]
real(DP), parameter :: pc     = 206265*AU                   !Astronomical unit         [cm]
real(DP), parameter :: pi=acos(-1.0d0), rho0 = 1e-4, x_star=0, y_star=0, r_inn = 1, r_out = 40 ![r] = au, [rho0] = cm**-3
real(DP), parameter :: pixel_size=5, d_obs = 150*pc ![pixel_size] = cm
real(DP) :: L_star, nu(3), Qsca(3), T_star = 5778, R_star = RSun
real(DP) :: r_g = 1.005E-05 ![r_g] = cm
real(DP), allocatable :: x(:), y(:)
real(DP), allocatable :: xb(:), yb(:), intensity(:,:), intensity1d(:)
real(DP), allocatable :: rho(:,:), a_sca(:,:), source_func(:,:)
integer :: i, j, inu, n_phot = 1e6, n_half_grid = 50, n_case
character(100) :: filename_source, filename_intens

!initialising grids
allocate(x(2*n_half_grid+1), y(2*n_half_grid+1))
allocate(xb(2*n_half_grid+2), yb(2*n_half_grid+2), intensity(2*n_half_grid+2, 2*n_half_grid+2), intensity1d(2*n_half_grid+1))
allocate(rho(2*n_half_grid+1, 2*n_half_grid+1), a_sca(2*n_half_grid+1, 2*n_half_grid+1), source_func(2*n_half_grid+1, 2*n_half_grid+1))

do i=1,2*n_half_grid+1
    x(i) = (i-(n_half_grid+1))
    y(i) = (i-(n_half_grid+1)) !in AU
enddo
do i=1,2*n_half_grid+2
    xb(i) = -0.5 + (i-(n_half_grid+1))
    yb(i) = -0.5 + (i-(n_half_grid+1)) !in AU
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n_case = 1 !1 - source function methon, 2 - monte carlo picture method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1,2*n_half_grid+1
    do i=1, 2*n_half_grid+1
        if ((x(i)**2+y(j)**2 .lt. r_out**2).and.(x(i)**2+y(j)**2 .gt. r_inn**2)) then
            rho(i,j) = rho0*r_inn**2/(x(i)**2+y(j)**2)
        else
            rho(i,j) = 0
        endif
    enddo
enddo
! a_sca(n_half_grid+1, n_half_grid+1) = 0
![nu] = s**-1, [Qsca] = dimensionless
nu(1) = clight/0.709E-04; Qsca(1) = 0.255E+00
nu(2) = clight/0.534E-04; Qsca(2) = 0.692E+00
nu(3) = clight/0.454E-04; Qsca(3) = 0.116E+01
Qsca = 0
do inu=1,3
    intensity = 0; intensity1d = 0
    source_func = 0
    a_sca = pi*r_g**2*Qsca(inu)*rho
    L_star = 2*pi*R_star*pi*B(nu(inu), T_star)
    select case(n_case)
    case(1)
        do i=1, n_phot
            call monte_carlo_sca()
            print*, i, ' photon finished'
        enddo
        write(filename_source, '(a, i1)') 'source_function0_', inu
        open(400, file=filename_source)
        do j=2*n_half_grid+1, 1, -1
            write(400, '(101(es12.3))') (source_func(i,j), i=1, 2*n_half_grid+1)
        enddo
        close(400)
        call integrate_transfer()
        write(filename_intens, '(a, i1)') 'intensity1d_raytrace0_', inu
        open(405, file=filename_intens)
        do j=1, 2*n_half_grid+1
            write(405, '(es12.3)') intensity(1, j)*AU/d_obs
        enddo
        close(405)
    case(2)
        do i=1, n_phot
            call monte_carlo_pic()
            print*, i, ' photon finished'
        enddo
        write(filename_intens, '(a, i1)') 'intensity1d_montecarlo_', inu
        open(203, file=filename_intens)
        do j=1, 2*n_half_grid+1
            write(203, '(es12.3)') intensity1d(j)
        enddo
        close(203)
    end select
enddo
! print '(es12.3)', minval(source_func)

contains

subroutine monte_carlo_sca()
implicit none
real(DP) :: x0, y0
real(DP) :: xi, tau, dtau, lum_pack, vol_c
real(DP) :: xi1, phi, cos_dpm
real(DP) :: dxm, dym, x1, y1, ds
integer :: im(1), jm(1), test

x0 = x_star; y0 = y_star
outer: do
!     print*, x0, y0
    call random_number(xi)
!     print*, xi
    tau = -log(xi)
!     print*, tau
    call random_number(xi1)
    phi = 2*pi*xi1
!     print*, phi*180/pi
    do
        if (minval(abs(xb-x0)).gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
        else if (cos(phi).gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).ge.0)
!             print*, 'a'
            if (x0.eq.xb(size(xb))) exit outer
        else
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
!             print*, x0, im, jm
            if (x0.eq.xb(1)) exit outer
        endif

        if (minval(abs(yb-y0)).gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
        else if (sin(phi).gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).ge.0)
            if (y0.eq.yb(size(yb))) exit outer
        else
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
            if (y0.eq.yb(1)) exit outer
        endif

!         print*, im, jm


        if (cos(phi).eq.0) then
            dxm = 0
            dym = sign(minval(abs(yb-y0), (yb-y0)*sin(phi).gt.0), sin(phi))
        else if (sin(phi).eq.0) then
            dxm = sign(minval(abs(xb-x0), (xb-x0)*cos(phi).gt.0), cos(phi))
            dym = 0
        else
            dxm = sign(minval(abs(xb-x0), (xb-x0)*cos(phi).gt.0), cos(phi))
            dym = sign(minval(abs(yb-y0), (yb-y0)*sin(phi).gt.0), sin(phi))
        endif

        cos_dpm = dxm/sqrt(dxm**2+dym**2)
        if (abs(cos(phi)).gt.abs(cos_dpm)) then
            x1 = x0 + dxm
            y1 = y0 + dxm*tan(phi)
        else if (abs(cos(phi)).lt.abs(cos_dpm)) then
            y1 = y0 + dym
            x1 = x0 + dym/tan(phi)
        else
            x1 = x0 + dxm
            y1 = y0 + dym
        endif
!         print*, dxm, dym
        ds = sqrt((x1-x0)**2 + (y1 - y0)**2)*AU ![ds] = cm
        dtau = a_sca(im(1), jm(1))*ds
        if (tau.le.dtau) ds = tau/a_sca(im(1), jm(1))
        tau = tau - dtau
        lum_pack = L_star/n_phot
        vol_c = abs(xb(im(1)+1)-xb(im(1)))*abs(yb(jm(1)+1)-yb(jm(1)))*AU**2 ![vol_c] = cm**2
        source_func(im(1),jm(1)) = source_func(im(1),jm(1)) + ds*lum_pack/4./pi/vol_c
!         print*, x0, y0
        if (tau.le.0) then
            x0 = x0 + ds*cos(phi)/AU
            y0 = y0 + ds*sin(phi)/AU
            exit
        else
            x0 = x1; y0 = y1
        endif
    enddo
enddo outer
end subroutine

subroutine monte_carlo_pic()
implicit none
real(DP) :: x0, y0
real(DP) :: xi, tau, dtau, lum_pack, vol_c
real(DP) :: xi1, phi, cos_dpm
real(DP) :: dxm, dym, x1, y1, ds
integer :: im(1), jm(1), test

x0 = x_star; y0 = y_star
outer: do
!     print*, x0, y0
    call random_number(xi)
!     print*, xi
    tau = -log(xi)
!     print*, tau
    call random_number(xi1)
    phi = 2*pi*xi1
!     print*, phi*180/pi
    do
        if (minval(abs(xb-x0)).gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
        else if (cos(phi).gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).ge.0)
!             print*, 'a'
            if (x0.eq.xb(size(xb))) exit outer
        else
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
!             print*, x0, im, jm
            if (x0.eq.xb(1)) then
                if ((cos(phi)+1.le.1e-5)) then
                   jm = minloc(abs(yb-y0), (yb-y0).lt.0)
                   lum_pack = L_star/n_phot
                   intensity1d(jm(1)) = intensity1d(jm(1)) + lum_pack/d_obs
                endif
                exit outer
            endif
        endif

        if (minval(abs(yb-y0)).gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
        else if (sin(phi).gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).ge.0)
            if (y0.eq.yb(size(yb))) exit outer
        else
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
            if (y0.eq.yb(1)) exit outer
        endif

!         print*, im, jm


        if (cos(phi).eq.0) then
            dxm = 0
            dym = sign(minval(abs(yb-y0), (yb-y0)*sin(phi).gt.0), sin(phi))
        else if (sin(phi).eq.0) then
            dxm = sign(minval(abs(xb-x0), (xb-x0)*cos(phi).gt.0), cos(phi))
            dym = 0
        else
            dxm = sign(minval(abs(xb-x0), (xb-x0)*cos(phi).gt.0), cos(phi))
            dym = sign(minval(abs(yb-y0), (yb-y0)*sin(phi).gt.0), sin(phi))
        endif

        cos_dpm = dxm/sqrt(dxm**2+dym**2)
        if (abs(cos(phi)).gt.abs(cos_dpm)) then
            x1 = x0 + dxm
            y1 = y0 + dxm*tan(phi)
        else if (abs(cos(phi)).lt.abs(cos_dpm)) then
            y1 = y0 + dym
            x1 = x0 + dym/tan(phi)
        else
            x1 = x0 + dxm
            y1 = y0 + dym
        endif
!         print*, dxm, dym
        ds = sqrt((x1-x0)**2 + (y1 - y0)**2)*AU
        dtau = a_sca(im(1), jm(1))*ds
        if (tau.le.dtau) ds = tau/a_sca(im(1), jm(1))
        tau = tau - dtau
        if (tau.le.0) then
            x0 = x0 + ds*cos(phi)/AU
            y0 = y0 + ds*sin(phi)/AU
            exit
        else
            x0 = x1; y0 = y1
        endif
    enddo
enddo outer
end subroutine

subroutine integrate_transfer()
implicit none
integer :: ii, jj
real(DP) :: dtau

do jj=1, 2*n_half_grid+1
    do ii=2*n_half_grid+1, 1, -1
        dtau =  a_sca(ii, jj)*(xb(ii+1)-xb(ii))*AU
        intensity(ii, jj) = intensity(ii+1, jj)*exp(-dtau) + (1 - exp(-dtau))*source_func(ii, jj)
        if ((jj.eq.n_half_grid+1).and.(ii.eq.n_half_grid+1)) intensity(ii,jj) = intensity(ii,jj) + L_star/(2*pi*R_star)*exp(-a_sca(ii,jj)*(x_star-xb(ii))*AU)
    enddo
enddo
end subroutine

function B(nu, T)
implicit none
real(DP) :: B, nu, T ![nu] = Hz, [T] = K, [B] = erg/(cm^2 s Hz)
B = 2*hP*nu**3/clight**2 * 1/(exp(hp*nu/kB/T)-1)
end function
end program part1
