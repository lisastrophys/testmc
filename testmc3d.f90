program part1_3d
!$ use omp_lib
implicit none
integer, PARAMETER :: DP = KIND(1.0d0)
real(DP), parameter :: RSun   = 6.957e10_DP                 !Solar radius              [cm]
real(DP), parameter :: clight = 2.99792458e10_DP            !Speed of light            [cm s-1]
real(DP), parameter :: hP     = 6.626070040e-27_DP          !Planck constant           [erg s]
real(DP), parameter :: kB     = 1.38064852e-16_DP           !Boltzmann constant        [erg K-1]
real(DP), parameter :: AU     = 1.495978707e13_DP           !Astronomical unit         [cm]
real(DP), parameter :: pc     = 206265*AU                   !Parsec                    [cm]
real(DP), parameter :: pi=acos(-1.0d0), rho0 = 1.2756e-18_DP, x_star=0, y_star=0, z_star=0, r_inn = 1, r_out = 40 ![r] = au, [rho0] = g cm**-3
real(DP), parameter :: pixel_size=5, d_obs = 150*pc ![pixel_size] = cm
real(DP), parameter :: dmu_obs = 1e-2   !half-width in cos(angle) of the cone towards the observer (-x) used in case 2
real(DP), parameter :: w_min = 1e-4, p_survive = 0.1 !russian roulette: packets with weight < w_min survive with probability p_survive
character(*), parameter :: dust_file = 'dustkappa.inp' !radmc3d dust opacity file
real(DP) :: T_star = 5778, R_star = RSun
real(DP), allocatable :: lam(:), nu(:), kappa_abs(:), kappa_sca(:), kappa_ext(:), albedo(:), L_star(:)
real(DP), allocatable :: x(:), y(:), z(:)
real(DP), allocatable :: xb(:), yb(:), zb(:)
real(DP), allocatable :: rho(:,:,:), source_func(:,:,:,:), image_rt(:,:,:), image_mc(:,:,:)
real(DP), allocatable :: temp(:,:,:) !dust temperature [K] (n_case = 3)
integer :: i, j, k, inu, nlam, n_phot = 1e6, n_half_grid = 50, n_case, n_cell
integer(8) :: ip, n_tot, n_report
character(100) :: filename
logical, allocatable :: out_lam(:) !.true. for wavelength indices whose results are written out

!reading dust optical properties
call read_dust_opac(dust_file)
!run mode and wavelength indices to write out from the command line, e.g. ./mc3d imagert inu 1 40 50
call parse_command_line()
allocate(nu(nlam), kappa_ext(nlam), albedo(nlam), L_star(nlam))
nu = clight/(lam*1e-4_DP)  ![lam] = micron, [nu] = s**-1
kappa_ext = kappa_abs + kappa_sca ![kappa] = cm**2 g**-1 (per gram of dust)
where (kappa_ext.gt.0)
    albedo = kappa_sca/kappa_ext
elsewhere
    albedo = 0
endwhere
do inu=1, nlam
    L_star(inu) = 4*pi*R_star**2*pi*B(nu(inu), T_star) ![L_star] = erg s-1 Hz-1
enddo

!initialising grids
n_cell = 2*n_half_grid+1
allocate(x(n_cell), y(n_cell), z(n_cell))
allocate(xb(n_cell+1), yb(n_cell+1), zb(n_cell+1))
allocate(rho(n_cell, n_cell, n_cell))
allocate(source_func(n_cell, n_cell, n_cell, nlam))
allocate(image_rt(n_cell, n_cell, nlam), image_mc(n_cell, n_cell, nlam))

do i=1,n_cell
    x(i) = (i-(n_half_grid+1))
    y(i) = (i-(n_half_grid+1))
    z(i) = (i-(n_half_grid+1)) !in AU
enddo
do i=1,n_cell+1
    xb(i) = -0.5 + (i-(n_half_grid+1))
    yb(i) = -0.5 + (i-(n_half_grid+1))
    zb(i) = -0.5 + (i-(n_half_grid+1)) !in AU
enddo
do k=1,n_cell
    do j=1,n_cell
        do i=1, n_cell
            if ((x(i)**2+y(j)**2+z(k)**2 .lt. r_out**2).and.(x(i)**2+y(j)**2+z(k)**2 .gt. r_inn**2)) then
                rho(i,j,k) = rho0*r_inn**2/(x(i)**2+y(j)**2+z(k)**2)
            else
                rho(i,j,k) = 0
            endif
        enddo
    enddo
enddo

!all wavelengths in one loop: packet ip has wavelength index mod(ip-1, nlam)+1, n_phot packets per wavelength
n_tot = int(n_phot,8)*nlam
n_report = max(n_tot/20, 1_8)
!$ print*, 'OpenMP threads: ', omp_get_max_threads()
print*, 'wavelengths: ', nlam, ' packets in total: ', n_tot
if (all(out_lam)) then
    print*, 'writing out wavelength indices: all'
else
    print*, 'writing out wavelength indices: ', pack([(inu, inu=1, nlam)], out_lam)
endif

select case(n_case)
case(1)
    source_func = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = int(mod(ip-1, int(nlam,8))) + 1
        call monte_carlo_sca(inu)
        if (mod(ip, n_report).eq.0) print*, ip, ' photons finished'
    enddo
    !$omp end parallel do
    !source_func holds the mean intensity J; scattering source function S = albedo*J
    do inu=1, nlam
        source_func(:,:,:,inu) = albedo(inu)*source_func(:,:,:,inu)
    enddo
    !full 3D cube per wavelength: blocks of constant k (k = 1..N), inside each block rows j = N..1, columns i = 1..N
    !do inu=1, nlam
    !    if (.not.out_lam(inu)) cycle
    !    write(filename, '(a, i4.4, a)') 'source_function3d_', inu, '.out'
    !    open(400, file=filename)
    !    do k=1, n_cell
    !        do j=n_cell, 1, -1
    !            write(400, '(*(es12.3))') (source_func(i,j,k,inu), i=1, n_cell)
    !        enddo
    !    enddo
    !    close(400)
    !enddo
    call integrate_transfer()
    !image seen from -x: rows k = N..1, columns j = 1..N, flux per pixel [erg s-1 cm-2 Hz-1]
    do inu=1, nlam
        if (.not.out_lam(inu)) cycle
        write(filename, '(a, i4.4, a)') 'image_raytrace_', inu, '.out'
        open(405, file=filename)
        do k=n_cell, 1, -1
            write(405, '(*(es12.3))') (image_rt(j, k, inu), j=1, n_cell)
        enddo
        close(405)
    enddo
case(2)
    image_mc = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = int(mod(ip-1, int(nlam,8))) + 1
        call monte_carlo_pic(inu)
        if (mod(ip, n_report).eq.0) print*, ip, ' photons finished'
    enddo
    !$omp end parallel do
    do inu=1, nlam
        if (.not.out_lam(inu)) cycle
        write(filename, '(a, i4.4, a)') 'image_montecarlo_', inu, '.out'
        open(203, file=filename)
        do k=n_cell, 1, -1
            write(203, '(*(es12.3))') (image_mc(j, k, inu), j=1, n_cell)
        enddo
        close(203)
    enddo
case(3)
    !mean intensity J from the same packets as in case 1 (no albedo factor: dust_temp needs J itself)
    source_func = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = int(mod(ip-1, int(nlam,8))) + 1
        call monte_carlo_sca(inu)
        if (mod(ip, n_report).eq.0) print*, ip, ' photons finished'
    enddo
    !$omp end parallel do
    allocate(temp(n_cell, n_cell, n_cell))
    call dust_temp()
    print*, 'dust temperature min/max [K]: ', minval(temp, rho.gt.0), maxval(temp, rho.gt.0)
    !full 3D cube: blocks of constant k (k = 1..N), inside each block rows j = N..1, columns i = 1..N
    open(410, file='dust_temperature3d.out')
    do k=1, n_cell
        do j=n_cell, 1, -1
            write(410, '(*(es12.3))') (temp(i,j,k), i=1, n_cell)
        enddo
    enddo
    close(410)
case(4)
    !mean intensity J from the same packets as in case 1 (no albedo factor: dust_temp needs J itself)
    source_func = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = int(mod(ip-1, int(nlam,8))) + 1
        call monte_carlo_sca(inu)
        if (mod(ip, n_report).eq.0) print*, ip, ' photons finished'
    enddo
    !$omp end parallel do
    allocate(temp(n_cell, n_cell, n_cell))
    call dust_temp()
    print*, 'dust temperature min/max [K]: ', minval(temp, rho.gt.0), maxval(temp, rho.gt.0)
    !full 3D cube: blocks of constant k (k = 1..N), inside each block rows j = N..1, columns i = 1..N
    open(410, file='dust_temperature3d.out')
    do k=1, n_cell
        do j=n_cell, 1, -1
            write(410, '(*(es12.3))') (temp(i,j,k), i=1, n_cell)
        enddo
    enddo
    close(410)
    do inu=1, nlam
        source_func(:,:,:,inu) = albedo(inu)*source_func(:,:,:,inu) + (1-albedo(inu))*B(nu(inu), temp)
    enddo
    call integrate_transfer()
    !image seen from -x: rows k = N..1, columns j = 1..N, flux per pixel [erg s-1 cm-2 Hz-1]
    do inu=1, nlam
        if (.not.out_lam(inu)) cycle
        write(filename, '(a, i4.4, a)') 'image_raytrace_therm_', inu, '.out'
        open(405, file=filename)
        do k=n_cell, 1, -1
            write(405, '(*(es12.3))') (image_rt(j, k, inu), j=1, n_cell)
        enddo
        close(405)
    enddo
end select

contains

subroutine monte_carlo_sca(inu)
implicit none
integer, intent(in) :: inu
real(DP) :: x0, y0, z0
real(DP) :: xi, tau, dtau, lum_pack, vol_c, w, a_ext
real(DP) :: xi1, xi2, phi, cos_theta, sin_theta, ux, uy, uz
real(DP) :: xw, yw, zw, sx, sy, sz, s_min, x1, y1, z1, ds
integer :: im(1), jm(1), km(1)

x0 = x_star; y0 = y_star; z0 = z_star
w = 1 !packet weight, reduced by the albedo at each scattering
lum_pack = L_star(inu)/n_phot
outer: do
    call random_number(xi)
    tau = -log(xi)
    !isotropic direction: phi uniform in [0, 2pi), cos(theta) uniform in [-1, 1]
    call random_number(xi1)
    call random_number(xi2)
    phi = 2*pi*xi1
    cos_theta = 2*xi2 - 1
    sin_theta = sqrt(1 - cos_theta**2)
    ux = sin_theta*cos(phi)
    uy = sin_theta*sin(phi)
    uz = cos_theta
    do
        !indices of the current cell (on a wall -> the cell the photon is moving into)
        if (minval(abs(xb-x0)).gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
        else if (ux.gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).ge.0)
            if (x0.eq.xb(size(xb))) exit outer
        else
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
            if (x0.eq.xb(1)) exit outer
        endif

        if (minval(abs(yb-y0)).gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
        else if (uy.gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).ge.0)
            if (y0.eq.yb(size(yb))) exit outer
        else
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
            if (y0.eq.yb(1)) exit outer
        endif

        if (minval(abs(zb-z0)).gt.0) then
            km = minloc(abs(zb-z0), (zb-z0).lt.0)
        else if (uz.gt.0) then
            km = minloc(abs(zb-z0), (zb-z0).ge.0)
            if (z0.eq.zb(size(zb))) exit outer
        else
            km = minloc(abs(zb-z0), (zb-z0).lt.0)
            if (z0.eq.zb(1)) exit outer
        endif

        !path length (in AU) to the x-, y- and z-walls of the current cell; the closest one is crossed
        if (ux.gt.0) then
            xw = xb(im(1)+1); sx = (xw-x0)/ux
        else if (ux.lt.0) then
            xw = xb(im(1));   sx = (xw-x0)/ux
        else
            xw = x0;          sx = huge(1.0_DP)
        endif
        if (uy.gt.0) then
            yw = yb(jm(1)+1); sy = (yw-y0)/uy
        else if (uy.lt.0) then
            yw = yb(jm(1));   sy = (yw-y0)/uy
        else
            yw = y0;          sy = huge(1.0_DP)
        endif
        if (uz.gt.0) then
            zw = zb(km(1)+1); sz = (zw-z0)/uz
        else if (uz.lt.0) then
            zw = zb(km(1));   sz = (zw-z0)/uz
        else
            zw = z0;          sz = huge(1.0_DP)
        endif

        s_min = min(sx, sy, sz)
        x1 = x0 + s_min*ux
        y1 = y0 + s_min*uy
        z1 = z0 + s_min*uz
        if (sx.eq.s_min) x1 = xw !put the photon exactly on the crossed wall
        if (sy.eq.s_min) y1 = yw
        if (sz.eq.s_min) z1 = zw

        ds = s_min*AU ![ds] = cm
        a_ext = rho(im(1), jm(1), km(1))*kappa_ext(inu) ![a_ext] = cm**-1
        dtau = a_ext*ds
        if (tau.le.dtau) ds = tau/a_ext
        tau = tau - dtau
        vol_c = abs(xb(im(1)+1)-xb(im(1)))*abs(yb(jm(1)+1)-yb(jm(1)))*abs(zb(km(1)+1)-zb(km(1)))*AU**3 ![vol_c] = cm**3
        !$omp atomic
        source_func(im(1),jm(1),km(1),inu) = source_func(im(1),jm(1),km(1),inu) + w*ds*lum_pack/4./pi/vol_c
        if (tau.le.0) then
            x0 = x0 + ds*ux/AU
            y0 = y0 + ds*uy/AU
            z0 = z0 + ds*uz/AU
            !interaction: absorbed fraction removed from the packet, the rest is scattered
            w = w*albedo(inu)
            if (w.lt.w_min) then
                call random_number(xi)
                if (xi.ge.p_survive) exit outer
                w = w/p_survive
            endif
            exit
        else
            x0 = x1; y0 = y1; z0 = z1
        endif
    enddo
enddo outer
end subroutine

subroutine monte_carlo_pic(inu)
implicit none
integer, intent(in) :: inu
real(DP) :: x0, y0, z0
real(DP) :: xi, tau, dtau, lum_pack, w, a_ext
real(DP) :: xi1, xi2, phi, cos_theta, sin_theta, ux, uy, uz
real(DP) :: xw, yw, zw, sx, sy, sz, s_min, x1, y1, z1, ds
real(DP) :: y_s, z_s !position of last emission/scattering, projected on the sky plane
integer :: im(1), jm(1), km(1)

x0 = x_star; y0 = y_star; z0 = z_star
w = 1
lum_pack = L_star(inu)/n_phot
outer: do
    call random_number(xi)
    tau = -log(xi)
    call random_number(xi1)
    call random_number(xi2)
    phi = 2*pi*xi1
    cos_theta = 2*xi2 - 1
    sin_theta = sqrt(1 - cos_theta**2)
    ux = sin_theta*cos(phi)
    uy = sin_theta*sin(phi)
    uz = cos_theta
    y_s = y0; z_s = z0
    do
        if (minval(abs(xb-x0)).gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
        else if (ux.gt.0) then
            im = minloc(abs(xb-x0), (xb-x0).ge.0)
            if (x0.eq.xb(size(xb))) exit outer
        else
            im = minloc(abs(xb-x0), (xb-x0).lt.0)
            if (x0.eq.xb(1)) then
                !photon leaves through the -x face inside the cone towards the observer
                if ((ux+1.le.dmu_obs)) then
                    jm = minloc(abs(yb-y_s), (yb-y_s).lt.0)
                    km = minloc(abs(zb-z_s), (zb-z_s).lt.0)
                    !energy per steradian of the cone (solid angle 2*pi*dmu_obs) divided by d_obs**2
                    if ((jm(1).ge.1).and.(km(1).ge.1)) then
                        !$omp atomic
                        image_mc(jm(1), km(1), inu) = image_mc(jm(1), km(1), inu) + w*lum_pack/(2*pi*dmu_obs)/d_obs**2
                    endif
                endif
                exit outer
            endif
        endif

        if (minval(abs(yb-y0)).gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
        else if (uy.gt.0) then
            jm = minloc(abs(yb-y0), (yb-y0).ge.0)
            if (y0.eq.yb(size(yb))) exit outer
        else
            jm = minloc(abs(yb-y0), (yb-y0).lt.0)
            if (y0.eq.yb(1)) exit outer
        endif

        if (minval(abs(zb-z0)).gt.0) then
            km = minloc(abs(zb-z0), (zb-z0).lt.0)
        else if (uz.gt.0) then
            km = minloc(abs(zb-z0), (zb-z0).ge.0)
            if (z0.eq.zb(size(zb))) exit outer
        else
            km = minloc(abs(zb-z0), (zb-z0).lt.0)
            if (z0.eq.zb(1)) exit outer
        endif

        if (ux.gt.0) then
            xw = xb(im(1)+1); sx = (xw-x0)/ux
        else if (ux.lt.0) then
            xw = xb(im(1));   sx = (xw-x0)/ux
        else
            xw = x0;          sx = huge(1.0_DP)
        endif
        if (uy.gt.0) then
            yw = yb(jm(1)+1); sy = (yw-y0)/uy
        else if (uy.lt.0) then
            yw = yb(jm(1));   sy = (yw-y0)/uy
        else
            yw = y0;          sy = huge(1.0_DP)
        endif
        if (uz.gt.0) then
            zw = zb(km(1)+1); sz = (zw-z0)/uz
        else if (uz.lt.0) then
            zw = zb(km(1));   sz = (zw-z0)/uz
        else
            zw = z0;          sz = huge(1.0_DP)
        endif

        s_min = min(sx, sy, sz)
        x1 = x0 + s_min*ux
        y1 = y0 + s_min*uy
        z1 = z0 + s_min*uz
        if (sx.eq.s_min) x1 = xw
        if (sy.eq.s_min) y1 = yw
        if (sz.eq.s_min) z1 = zw

        ds = s_min*AU
        a_ext = rho(im(1), jm(1), km(1))*kappa_ext(inu)
        dtau = a_ext*ds
        if (tau.le.dtau) ds = tau/a_ext
        tau = tau - dtau
        if (tau.le.0) then
            x0 = x0 + ds*ux/AU
            y0 = y0 + ds*uy/AU
            z0 = z0 + ds*uz/AU
            w = w*albedo(inu)
            if (w.lt.w_min) then
                call random_number(xi)
                if (xi.ge.p_survive) exit outer
                w = w/p_survive
            endif
            exit
        else
            x0 = x1; y0 = y1; z0 = z1
        endif
    enddo
enddo outer
end subroutine

subroutine integrate_transfer()
implicit none
integer :: ii, jj, kk, inu
real(DP) :: dtau, a_ext, I_ray

!rays along -x towards the observer, one per (y, z) pixel and wavelength
!$omp parallel do collapse(3) private(ii, dtau, a_ext, I_ray) schedule(static)
do inu=1, nlam
    do kk=1, n_cell
        do jj=1, n_cell
            I_ray = 0
            do ii=n_cell, 1, -1
                a_ext = rho(ii, jj, kk)*kappa_ext(inu)
                dtau = a_ext*(xb(ii+1)-xb(ii))*AU
                I_ray = I_ray*exp(-dtau) + (1 - exp(-dtau))*source_func(ii, jj, kk, inu)
                !unresolved star: its flux spread over the pixel area -> mean intensity of the central pixel
                if ((kk.eq.n_half_grid+1).and.(jj.eq.n_half_grid+1).and.(ii.eq.n_half_grid+1)) &
                    I_ray = I_ray + L_star(inu)/(4*pi*(yb(jj+1)-yb(jj))*(zb(kk+1)-zb(kk))*AU**2) &
                                    *exp(-a_ext*(x_star-xb(ii))*AU)
            enddo
            image_rt(jj, kk, inu) = I_ray*(AU/d_obs)**2 !flux per pixel
        enddo
    enddo
enddo
!$omp end parallel do
end subroutine

subroutine parse_command_line()
!usage: ./mc3d [imagert | imagemc | mctherm] [inu i1 i2 ...] [nphot N]
!   imagert - n_case = 1, source function method + ray-traced images (default)
!   imagemc - n_case = 2, monte carlo picture method
!   mctherm - n_case = 3, dust temperature
!   inu     - wavelength indices to write out (default: all)
!   nphot   - number of photon packets per wavelength, e.g. 100000 or 1e5 (default: value in the code)
implicit none
character(100) :: arg
integer :: iarg, nargs, ival, ios
real(DP) :: rval
logical :: reading_inu, inu_given, reading_nphot

allocate(out_lam(nlam))
out_lam = .false.
n_case = 0
reading_inu = .false.; inu_given = .false.; reading_nphot = .false.
nargs = command_argument_count()
do iarg=1, nargs
    call get_command_argument(iarg, arg)
    select case(trim(arg))
    case('imagert', 'imagemc', 'mctherm', 'imagetherm')
        if (n_case.ne.0) call usage('only one of imagert, imagemc, mctherm, imagetherm can be given')
        if (trim(arg).eq.'imagert') n_case = 1
        if (trim(arg).eq.'imagemc') n_case = 2
        if (trim(arg).eq.'mctherm') n_case = 3
        if (trim(arg).eq.'imagetherm') n_case = 4
        reading_inu = .false.
        if (reading_nphot) call usage('no number given after "nphot"')
    case('inu')
        if (reading_nphot) call usage('no number given after "nphot"')
        reading_inu = .true.; inu_given = .true.
    case('nphot')
        if (reading_nphot) call usage('no number given after "nphot"')
        reading_inu = .false.; reading_nphot = .true.
    case('-h', '--help', 'help')
        call usage('')
    case default
        if (reading_nphot) then
            !read as real so that 1e6 is accepted as well as 1000000
            read(arg, *, iostat=ios) rval
            if (ios.ne.0) call usage('not a number of photon packets: '//trim(arg))
            if ((rval.lt.1).or.(rval.gt.huge(n_phot)).or.(rval.ne.anint(rval))) &
                call usage('nphot must be a whole number between 1 and 2147483647: '//trim(arg))
            n_phot = nint(rval)
            reading_nphot = .false.
            cycle
        endif
        if (.not.reading_inu) call usage('unknown argument: '//trim(arg))
        read(arg, *, iostat=ios) ival
        if (ios.ne.0) call usage('not an integer wavelength index: '//trim(arg))
        if ((ival.lt.1).or.(ival.gt.nlam)) then
            print*, 'wavelength index ', ival, ' out of range 1 ..', nlam
            stop
        endif
        out_lam(ival) = .true.
    end select
enddo
if (reading_nphot) call usage('no number given after "nphot"')
if (n_case.eq.0) then
    n_case = 1
    print*, 'no run mode given, using imagert'
endif
if (.not.inu_given) then
    out_lam = .true.
else if (.not.any(out_lam)) then
    call usage('no wavelength indices given after "inu"')
endif
select case(n_case)
case(1); print*, 'run mode: imagert (source function method + ray tracing)'
case(2); print*, 'run mode: imagemc (monte carlo picture method)'
case(3); print*, 'run mode: mctherm (dust temperature)'
end select
print*, 'photon packets per wavelength: ', n_phot
end subroutine

subroutine usage(msg)
implicit none
character(*), intent(in) :: msg
if (len_trim(msg).gt.0) print*, trim(msg)
print*, 'usage: ./mc3d [imagert | imagemc | mctherm] [inu i1 i2 ...] [nphot N]'
print*, '   imagert - source function method + ray-traced images (default)'
print*, '   imagemc - monte carlo picture method'
print*, '   mctherm - dust temperature'
print*, '   inu     - wavelength indices to write out (default: all)'
print*, '   nphot   - photon packets per wavelength, e.g. 100000 or 1e5'
stop
end subroutine

subroutine read_dust_opac(fname)
!radmc3d dustkappa_*.inp format, lines starting with '#' or '!' are comments:
!   iformat     (1: lambda, kappa_abs   2: + kappa_sca   3: + g)
!   nlam
!   lambda[micron]  kappa_abs[cm2/g]  (kappa_sca[cm2/g])  (g)    nlam lines
!g is read but not used: scattering is treated as isotropic
implicit none
character(*), intent(in) :: fname
character(512) :: line
integer :: iformat, il

open(10, file=fname, status='old', action='read')
call read_data_line(10, line); read(line, *) iformat
if ((iformat.lt.1).or.(iformat.gt.3)) stop 'read_dust_opac: iformat must be 1, 2 or 3'
call read_data_line(10, line); read(line, *) nlam
allocate(lam(nlam), kappa_abs(nlam), kappa_sca(nlam))
kappa_sca = 0
do il=1, nlam
    call read_data_line(10, line)
    if (iformat.eq.1) then
        read(line, *) lam(il), kappa_abs(il)
    else
        read(line, *) lam(il), kappa_abs(il), kappa_sca(il)
    endif
enddo
close(10)
end subroutine

subroutine read_data_line(unit, line)
implicit none
integer, intent(in) :: unit
character(*), intent(out) :: line
integer :: ios
do
    read(unit, '(a)', iostat=ios) line
    if (ios.ne.0) stop 'read_dust_opac: unexpected end of file'
    line = adjustl(line)
    if (len_trim(line).eq.0) cycle
    if ((line(1:1).eq.'#').or.(line(1:1).eq.'!')) cycle
    exit
enddo
end subroutine

subroutine dust_temp()
!dust temperature from radiative equilibrium in every cell (Lucy 1999):
!   int kappa_abs J_nu dnu = int kappa_abs B_nu(T) dnu
!J_nu is taken from source_func filled by monte_carlo_sca (direct + scattered starlight).
!Heating by the dust's own thermal emission is not included (dust optically thin in the IR).
!The frequency integrals use trapezoidal weights on the wavelength grid of the opacity file,
!which must be monotonic and cover both the stellar spectrum and the dust emission.
implicit none
integer, parameter :: n_T = 2000
real(DP), parameter :: T_min = 1, T_max = 1e5 ![T] = K, range of the emission table
real(DP) :: T_tab(n_T), E_tab(n_T), dnu(nlam), E_abs
integer :: il, it, ii, jj, kk, lo, hi, mid

temp = 0
if (maxval(kappa_abs).le.0) then
    print*, 'dust_temp: kappa_abs = 0 at all wavelengths, no dust heating, temperature set to 0'
    return
endif
if (nlam.lt.2) stop 'dust_temp: at least 2 wavelengths are needed'
if (.not.(all(nu(2:nlam).gt.nu(1:nlam-1)).or.all(nu(2:nlam).lt.nu(1:nlam-1)))) &
    stop 'dust_temp: wavelength grid must be monotonic'

!trapezoidal frequency weights
dnu(1)    = 0.5*abs(nu(2)-nu(1))
dnu(nlam) = 0.5*abs(nu(nlam)-nu(nlam-1))
do il=2, nlam-1
    dnu(il) = 0.5*abs(nu(il+1)-nu(il-1))
enddo

!emission table E(T) = sum kappa_abs B(T) dnu, log-spaced in T
do it=1, n_T
    T_tab(it) = T_min*(T_max/T_min)**(real(it-1,DP)/(n_T-1))
    E_tab(it) = sum(kappa_abs*B(nu, T_tab(it))*dnu)
enddo

!$omp parallel do collapse(3) private(E_abs, lo, hi, mid) schedule(static)
do kk=1, n_cell
    do jj=1, n_cell
        do ii=1, n_cell
            if (rho(ii,jj,kk).le.0) cycle !no dust -> temperature left at 0
            E_abs = sum(kappa_abs*source_func(ii,jj,kk,:)*dnu)
            if (E_abs.le.0) then
                temp(ii,jj,kk) = 0
            else if (E_abs.le.E_tab(1)) then
                temp(ii,jj,kk) = T_min
            else if (E_abs.ge.E_tab(n_T)) then
                temp(ii,jj,kk) = T_max
            else
                !bisection in the table, then linear interpolation
                lo = 1; hi = n_T
                do while (hi-lo.gt.1)
                    mid = (lo+hi)/2
                    if (E_tab(mid).le.E_abs) then
                        lo = mid
                    else
                        hi = mid
                    endif
                enddo
                temp(ii,jj,kk) = T_tab(lo) + (E_abs-E_tab(lo))/(E_tab(hi)-E_tab(lo))*(T_tab(hi)-T_tab(lo))
            endif
        enddo
    enddo
enddo
!$omp end parallel do
end subroutine

elemental function B(nu, T)
implicit none
real(DP), intent(in) :: nu, T ![nu] = Hz, [T] = K, [B] = erg/(cm^2 s Hz sr)
real(DP) :: B, xx
xx = hP*nu/kB/T
if (xx.gt.700) then !avoid overflow of exp, B is negligible there
    B = 0
else
    B = 2*hP*nu**3/clight**2 * 1/(exp(xx)-1)
endif
end function
end program part1_3d
