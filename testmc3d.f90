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
real(DP), parameter :: LSun   = 3.828e33_DP                 !Solar luminosity          [erg s-1]
real(DP), parameter :: sigma_SB = 5.670374419e-5_DP         !Stefan-Boltzmann constant [erg cm-2 s-1 K-4]
real(DP), parameter :: T_cmb  = 2.73_DP                     !CMB temperature, background dust heating [K]
real(DP), parameter :: pi=acos(-1.0d0), x_star=0, y_star=0, z_star=0 ![x_star] = AU
!default model (no grid.inp/dust_density.inp): spherical shell rho = rho0*(r_inn/r)**2 for r_inn < r < r_out
real(DP), parameter :: rho0 = 1.2756e-18_DP, r_inn = 1, r_out = 40 ![r] = AU, [rho0] = g cm**-3
character(*), parameter :: grid_file = 'grid.inp', dens_file = 'dust_density.inp'
real(DP), parameter :: d_obs = 150*pc !distance of the observer [cm]
real(DP), parameter :: dmu_obs = 1e-2   !half-width in cos(angle) of the cone towards the observer (-x) used in case 2
!modified random walk (Min et al. 2009): used in cells where the sphere to the nearest wall has
!optical depth > mrw_gamma; the scatterings inside the sphere are replaced by one diffusion step
real(DP), parameter :: mrw_gamma = 10
real(DP), parameter :: w_min = 1e-4, p_survive = 0.1 !russian roulette: packets with weight < w_min survive with probability p_survive
character(*), parameter :: dust_file = 'dustkappa.inp' !radmc3d dust opacity file
character(*), parameter :: lc_file = 'lightcurve.inp'  !stellar light curve; if present the run is time-dependent
real(DP), parameter :: n_burn = 1.5 !equilibrium start: packet ages up to n_burn light crossing times of the grid diagonal
integer, parameter :: f_eq_max = 20 !equilibrium start: at most f_eq_max*nphot packets per wavelength
!photon packet in flight (time-dependent run)
type packet
    real(DP) :: x, y, z, ux, uy, uz !position [AU] and direction
    real(DP) :: tau                 !optical depth left to the next interaction
    real(DP) :: w, e                !weight and energy [erg Hz-1] (energy deposited = e*w)
    real(DP) :: t                   !time [s]
    real(DP) :: ys, zs, tobs        !sky position and observer time of the last emission/scattering
    real(DP) :: l_walk, r_walk      !modified random walk in progress: remaining path length [cm], sphere radius [AU]
    integer :: inu
    integer :: ix, iy, iz           !current cell (0 or n+1: left the grid)
    logical :: heat                 !heating packet (wavelength drawn from the stellar spectrum)
end type packet
real(DP) :: T_star = 5778, R_star = RSun
real(DP), allocatable :: lam(:), nu(:), kappa_abs(:), kappa_sca(:), kappa_ext(:), albedo(:), L_star(:)
real(DP), allocatable :: x(:), y(:), z(:)
real(DP), allocatable :: xb(:), yb(:), zb(:)
real(DP), allocatable :: rho(:,:,:), source_func(:,:,:,:), image_rt(:,:,:), image_mc(:,:,:)
real(DP), allocatable :: temp(:,:,:) !dust temperature [K] (n_case = 3)
integer :: i, j, k, inu, nlam, n_phot = 1e6, n_half_grid = 50, n_case
integer :: nx, ny, nz                    !number of cells along x, y, z
integer :: is_star, js_star, ks_star     !cell containing the star
integer(8) :: ip, n_tot, n_report
character(100) :: filename
logical, allocatable :: out_lam(:) !.true. for wavelength indices whose results are written out
!time-dependent run
logical :: timedep, eq_start, thermal, images_rt, images_mc
integer :: lc_iformat, lc_n, r_first, n_sel, n_act
integer, allocatable :: lc_it(:), sel(:), isel_of(:), act(:)
real(DP), allocatable :: lc_t(:), lc_Lb(:), lc_Rs(:), lc_Ts(:), tbin_lo(:), tbin_hi(:)
real(DP), allocatable :: w_rt(:,:,:,:), star_att(:), delta_x(:), img_t(:,:,:,:)
real(DP), allocatable :: s_prev(:,:,:,:), s_cur(:,:,:,:) !source function of the previous and current step (selected wavelengths)
real(DP) :: t_init, t_burn, inv_dt, c_prev !c_prev: time the previous step average is assigned to
type(packet), allocatable :: pk_old(:), pk_new(:) !census of packets in flight
real(DP), allocatable :: y_mrw(:), p_mrw(:) !modified random walk: probability p_mrw(y) to be still inside the sphere
integer(8) :: n_new, n_old_census = 0
!heating packets (mctherm, imagetherm): wavelengths drawn from the stellar spectrum, all with the same energy
real(DP), allocatable :: heat_abs(:,:,:)            !absorbed power per cell: sum kappa_abs J_nu dnu
real(DP), allocatable :: dnu_w(:), cdf_heat(:), p_heat(:) !frequency weights, cumulative and single probabilities
real(DP) :: L_heat                                  !stellar luminosity integrated over the wavelength grid [erg s-1]
logical :: inu_given                                !inu given on the command line
real(DP) :: mem_avail, mem_static                   !available memory (0: unknown) and memory of the arrays [bytes]
integer(8) :: n_long_walk !equilibrium run: random walks lasting longer than t_burn

!help works without input files
do i=1, command_argument_count()
    call get_command_argument(i, filename)
    if ((trim(filename).eq.'-h').or.(trim(filename).eq.'--help').or.(trim(filename).eq.'help')) call usage('')
enddo
!reading dust optical properties
call read_dust_opac(dust_file)
!run mode and wavelength indices to write out from the command line, e.g. ./mc3d imagert inu 1 40 50
call parse_command_line()
!selected wavelengths: images and J dumps, each gets nphot packets
n_sel = count(out_lam)
allocate(sel(n_sel), isel_of(nlam))
sel = pack([(i, i=1, nlam)], out_lam)
isel_of = 0
do i=1, n_sel
    isel_of(sel(i)) = i
enddo
thermal = (n_case.eq.3).or.(n_case.eq.4)
images_rt = (n_case.eq.1).or.(n_case.eq.4)
images_mc = (n_case.eq.2)
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
call mrw_init()
!time-dependent run if lightcurve.inp is in the folder
inquire(file=lc_file, exist=timedep)
if (timedep) call read_lightcurve(lc_file)

!grid and dust density: from grid.inp and dust_density.inp, or the default spherical shell
call setup_grid_density()
allocate(source_func(nx, ny, nz, n_sel)) !J (then S) of the selected wavelengths only
if (thermal) then
    allocate(heat_abs(nx, ny, nz), temp(nx, ny, nz))
    call setup_heating()
endif
call memory_check()
allocate(image_rt(ny, nz, nlam), image_mc(ny, nz, nlam))

!all wavelengths in one loop: packet ip has wavelength index mod(ip-1, nlam)+1, n_phot packets per wavelength
n_tot = int(n_phot,8)*nlam
n_report = max(n_tot/20, 1_8)
!$ print*, 'OpenMP threads: ', omp_get_max_threads()
print*, 'wavelengths: ', nlam
if (.not.any(out_lam)) then
    print*, 'writing out wavelength indices: none'
else if (all(out_lam)) then
    print*, 'writing out wavelength indices: all'
else
    print*, 'writing out wavelength indices: ', pack([(inu, inu=1, nlam)], out_lam)
endif

if (timedep) then
    call run_timedep()
else
!packets: nphot per selected wavelength, plus nphot heating packets in the thermal modes
n_act = n_sel
allocate(act(n_act)); act = sel
n_tot = int(n_phot,8)*n_act
n_report = max(n_tot/20, 1_8)
print*, 'packets: ', n_tot, ' for the selected wavelengths', merge(int(n_phot,8), 0_8, thermal), ' for heating'
select case(n_case)
case(1)
    source_func = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = act(int(mod(ip-1, int(n_act,8))) + 1)
        call monte_carlo_sca(inu, .false.)
        if (mod(ip, n_report).eq.0) print*, ip, ' photons finished'
    enddo
    !$omp end parallel do
    !source_func holds the mean intensity J; scattering source function S = albedo*J
    do i=1, n_sel
        source_func(:,:,:,i) = albedo(sel(i))*source_func(:,:,:,i)
    enddo
    !full 3D cube per wavelength: blocks of constant k (k = 1..N), inside each block rows j = N..1, columns i = 1..N
    !do inu=1, nlam
    !    if (.not.out_lam(inu)) cycle
    !    write(filename, '(a, i4.4, a)') 'source_function3d_', inu, '.out'
    !    open(400, file=filename)
    !    do k=1, nz
    !        do j=ny, 1, -1
    !            write(400, '(*(es13.3e3))') (source_func(i,j,k,inu), i=1, nx)
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
        do k=nz, 1, -1
            write(405, '(*(es13.3e3))') (image_rt(j, k, inu), j=1, ny)
        enddo
        close(405)
    enddo
case(2)
    image_mc = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = act(int(mod(ip-1, int(n_act,8))) + 1)
        call monte_carlo_pic(inu)
        if (mod(ip, n_report).eq.0) print*, ip, ' photons finished'
    enddo
    !$omp end parallel do
    do inu=1, nlam
        if (.not.out_lam(inu)) cycle
        write(filename, '(a, i4.4, a)') 'image_montecarlo_', inu, '.out'
        open(203, file=filename)
        do k=nz, 1, -1
            write(203, '(*(es13.3e3))') (image_mc(j, k, inu), j=1, ny)
        enddo
        close(203)
    enddo
case(3)
    !heating packets: wavelengths drawn from the stellar spectrum
    call set_heat_cdf(L_star)
    heat_abs = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_phot
        inu = draw_heat_lambda()
        call monte_carlo_sca(inu, .true.)
        call progress('heating packets', ip, int(n_phot,8))
    enddo
    !$omp end parallel do
    call dust_temp()
    print*, 'dust temperature min/max [K]: ', minval(temp, rho.gt.0), maxval(temp, rho.gt.0)
    !full 3D cube: blocks of constant k (k = 1..N), inside each block rows j = N..1, columns i = 1..N
    open(410, file='dust_temperature3d.out')
    do k=1, nz
        do j=ny, 1, -1
            write(410, '(*(es13.3e3))') (temp(i,j,k), i=1, nx)
        enddo
    enddo
    close(410)
case(4)
    !heating packets: wavelengths drawn from the stellar spectrum
    call set_heat_cdf(L_star)
    heat_abs = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_phot
        inu = draw_heat_lambda()
        call monte_carlo_sca(inu, .true.)
        call progress('heating packets', ip, int(n_phot,8))
    enddo
    !$omp end parallel do
    !J of the selected wavelengths for the scattered light in the images
    source_func = 0
    !$omp parallel do schedule(dynamic, 100) private(inu)
    do ip=1, n_tot
        inu = act(int(mod(ip-1, int(n_act,8))) + 1)
        call monte_carlo_sca(inu, .false.)
        call progress('image wavelength packets', ip, n_tot)
    enddo
    !$omp end parallel do
    call dust_temp()
    print*, 'dust temperature min/max [K]: ', minval(temp, rho.gt.0), maxval(temp, rho.gt.0)
    !full 3D cube: blocks of constant k (k = 1..N), inside each block rows j = N..1, columns i = 1..N
    open(410, file='dust_temperature3d.out')
    do k=1, nz
        do j=ny, 1, -1
            write(410, '(*(es13.3e3))') (temp(i,j,k), i=1, nx)
        enddo
    enddo
    close(410)
    do i=1, n_sel
        source_func(:,:,:,i) = albedo(sel(i))*source_func(:,:,:,i) + (1-albedo(sel(i)))*B(nu(sel(i)), temp)
    enddo
    call integrate_transfer()
    !image seen from -x: rows k = N..1, columns j = 1..N, flux per pixel [erg s-1 cm-2 Hz-1]
    do inu=1, nlam
        if (.not.out_lam(inu)) cycle
        write(filename, '(a, i4.4, a)') 'image_raytrace_therm_', inu, '.out'
        open(405, file=filename)
        do k=nz, 1, -1
            write(405, '(*(es13.3e3))') (image_rt(j, k, inu), j=1, ny)
        enddo
        close(405)
    enddo
end select
endif

contains

subroutine monte_carlo_sca(inu, heating)
!heating = .true.: heating packet (energy L_heat/nphot), adds kappa_abs*J to heat_abs
!heating = .false.: packet of a selected wavelength (energy L_nu/nphot), adds J to source_func
implicit none
integer, intent(in) :: inu
logical, intent(in) :: heating
real(DP) :: x0, y0, z0
real(DP) :: xi, tau, dtau, lum_pack, vol_c, w, a_ext
real(DP) :: xi1, xi2, phi, cos_theta, sin_theta, ux, uy, uz
real(DP) :: xw, yw, zw, sx, sy, sz, s_min, x1, y1, z1, ds
real(DP) :: r0, l_w, a_abs, ex, ey, ez !modified random walk
integer :: im(1), jm(1), km(1)
logical :: new_flight

x0 = x_star; y0 = y_star; z0 = z_star
w = 1 !packet weight, reduced by the albedo at each scattering
if (heating) then
    lum_pack = L_heat/n_phot
else
    lum_pack = L_star(inu)/n_phot
endif
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
    new_flight = .true.
    do
        !current cell: binary search at the start of a flight, afterwards updated at every wall crossing
        if (new_flight) then
            im(1) = cell_index(xb, x0, ux); jm(1) = cell_index(yb, y0, uy); km(1) = cell_index(zb, z0, uz)
            new_flight = .false.
        endif
        if ((im(1).lt.1).or.(im(1).gt.nx).or.(jm(1).lt.1).or.(jm(1).gt.ny).or.(km(1).lt.1).or.(km(1).gt.nz)) exit outer

        !modified random walk in optically thick cells
        a_ext = rho(im(1), jm(1), km(1))*kappa_ext(inu)
        r0 = min(x0-xb(im(1)), xb(im(1)+1)-x0, y0-yb(jm(1)), yb(jm(1)+1)-y0, z0-zb(km(1)), zb(km(1)+1)-z0) ![r0] = AU
        if (a_ext*r0*AU.gt.mrw_gamma) then
            l_w = mrw_length(a_ext, r0*AU) ![l_w] = cm
            a_abs = rho(im(1), jm(1), km(1))*kappa_abs(inu)
            vol_c = abs(xb(im(1)+1)-xb(im(1)))*abs(yb(jm(1)+1)-yb(jm(1)))*abs(zb(km(1)+1)-zb(km(1)))*AU**3
            call deposit(im(1), jm(1), km(1), inu, heating, w*mrw_eff_length(a_abs, l_w)*lum_pack/4./pi/vol_c)
            w = w*exp(-a_abs*l_w)
            call random_dir(ex, ey, ez)
            x0 = x0 + r0*ex; y0 = y0 + r0*ey; z0 = z0 + r0*ez
            !the sphere lies inside the cell: keep rounding errors from moving the photon out of it
            x0 = min(max(x0, xb(im(1))), xb(im(1)+1))
            y0 = min(max(y0, yb(jm(1))), yb(jm(1)+1))
            z0 = min(max(z0, zb(km(1))), zb(km(1)+1))
            if (w.lt.w_min) then
                call random_number(xi)
                if (xi.ge.p_survive) exit outer
                w = w/p_survive
            endif
            !leaves the sphere surface outwards (cosine-weighted around the normal), new optical depth
            call outward_dir(ex, ey, ez, ux, uy, uz)
            call random_number(xi)
            tau = -log(xi)
            cycle
        endif

        !path length (in AU) to the x-, y- and z-walls of the current cell; the closest one is crossed
        if (ux.gt.0) then
            xw = xb(im(1)+1); sx = max((xw-x0)/ux, 0.0_DP)
        else if (ux.lt.0) then
            xw = xb(im(1));   sx = max((xw-x0)/ux, 0.0_DP)
        else
            xw = x0;          sx = huge(1.0_DP)
        endif
        if (uy.gt.0) then
            yw = yb(jm(1)+1); sy = max((yw-y0)/uy, 0.0_DP)
        else if (uy.lt.0) then
            yw = yb(jm(1));   sy = max((yw-y0)/uy, 0.0_DP)
        else
            yw = y0;          sy = huge(1.0_DP)
        endif
        if (uz.gt.0) then
            zw = zb(km(1)+1); sz = max((zw-z0)/uz, 0.0_DP)
        else if (uz.lt.0) then
            zw = zb(km(1));   sz = max((zw-z0)/uz, 0.0_DP)
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
        call deposit(im(1), jm(1), km(1), inu, heating, w*ds*lum_pack/4./pi/vol_c)
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
            !wall crossing: neighbouring cell
            if (sx.eq.s_min) im(1) = im(1) + merge(1, -1, ux.gt.0)
            if (sy.eq.s_min) jm(1) = jm(1) + merge(1, -1, uy.gt.0)
            if (sz.eq.s_min) km(1) = km(1) + merge(1, -1, uz.gt.0)
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
real(DP) :: r0, l_w, a_abs, ex, ey, ez !modified random walk
real(DP) :: y_s, z_s !position of last emission/scattering, projected on the sky plane
integer :: im(1), jm(1), km(1)
logical :: new_flight

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
    new_flight = .true.
    do
        !current cell: binary search at the start of a flight, afterwards updated at every wall crossing
        if (new_flight) then
            im(1) = cell_index(xb, x0, ux); jm(1) = cell_index(yb, y0, uy); km(1) = cell_index(zb, z0, uz)
            new_flight = .false.
        endif
        if (im(1).lt.1) then
            !photon leaves through the -x face inside the cone towards the observer
            if ((ux+1.le.dmu_obs)) then
                jm(1) = cell_index(yb, y_s, -1.0_DP)
                km(1) = cell_index(zb, z_s, -1.0_DP)
                !energy per steradian of the cone (solid angle 2*pi*dmu_obs) divided by d_obs**2
                if ((jm(1).ge.1).and.(jm(1).le.ny).and.(km(1).ge.1).and.(km(1).le.nz)) then
                    !$omp atomic
                    image_mc(jm(1), km(1), inu) = image_mc(jm(1), km(1), inu) + w*lum_pack/(2*pi*dmu_obs)/d_obs**2
                endif
            endif
            exit outer
        endif
        if ((im(1).gt.nx).or.(jm(1).lt.1).or.(jm(1).gt.ny).or.(km(1).lt.1).or.(km(1).gt.nz)) exit outer

        !modified random walk in optically thick cells
        a_ext = rho(im(1), jm(1), km(1))*kappa_ext(inu)
        r0 = min(x0-xb(im(1)), xb(im(1)+1)-x0, y0-yb(jm(1)), yb(jm(1)+1)-y0, z0-zb(km(1)), zb(km(1)+1)-z0)
        if (a_ext*r0*AU.gt.mrw_gamma) then
            l_w = mrw_length(a_ext, r0*AU)
            a_abs = rho(im(1), jm(1), km(1))*kappa_abs(inu)
            w = w*exp(-a_abs*l_w)
            call random_dir(ex, ey, ez)
            x0 = x0 + r0*ex; y0 = y0 + r0*ey; z0 = z0 + r0*ez
            !the sphere lies inside the cell: keep rounding errors from moving the photon out of it
            x0 = min(max(x0, xb(im(1))), xb(im(1)+1))
            y0 = min(max(y0, yb(jm(1))), yb(jm(1)+1))
            z0 = min(max(z0, zb(km(1))), zb(km(1)+1))
            if (w.lt.w_min) then
                call random_number(xi)
                if (xi.ge.p_survive) exit outer
                w = w/p_survive
            endif
            call outward_dir(ex, ey, ez, ux, uy, uz)
            call random_number(xi)
            tau = -log(xi)
            y_s = y0; z_s = z0 !the walk acts like the last scattering
            cycle
        endif

        if (ux.gt.0) then
            xw = xb(im(1)+1); sx = max((xw-x0)/ux, 0.0_DP)
        else if (ux.lt.0) then
            xw = xb(im(1));   sx = max((xw-x0)/ux, 0.0_DP)
        else
            xw = x0;          sx = huge(1.0_DP)
        endif
        if (uy.gt.0) then
            yw = yb(jm(1)+1); sy = max((yw-y0)/uy, 0.0_DP)
        else if (uy.lt.0) then
            yw = yb(jm(1));   sy = max((yw-y0)/uy, 0.0_DP)
        else
            yw = y0;          sy = huge(1.0_DP)
        endif
        if (uz.gt.0) then
            zw = zb(km(1)+1); sz = max((zw-z0)/uz, 0.0_DP)
        else if (uz.lt.0) then
            zw = zb(km(1));   sz = max((zw-z0)/uz, 0.0_DP)
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
            !wall crossing: neighbouring cell
            if (sx.eq.s_min) im(1) = im(1) + merge(1, -1, ux.gt.0)
            if (sy.eq.s_min) jm(1) = jm(1) + merge(1, -1, uy.gt.0)
            if (sz.eq.s_min) km(1) = km(1) + merge(1, -1, uz.gt.0)
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
    do kk=1, nz
        do jj=1, ny
            if (.not.out_lam(inu)) cycle
            I_ray = 0
            do ii=nx, 1, -1
                a_ext = rho(ii, jj, kk)*kappa_ext(inu)
                dtau = a_ext*(xb(ii+1)-xb(ii))*AU
                I_ray = I_ray*exp(-dtau) + (1 - exp(-dtau))*source_func(ii, jj, kk, isel_of(inu))
                !unresolved star: its flux spread over the pixel area -> mean intensity of the central pixel
                if ((kk.eq.ks_star).and.(jj.eq.js_star).and.(ii.eq.is_star)) &
                    I_ray = I_ray + L_star(inu)/(4*pi*(yb(jj+1)-yb(jj))*(zb(kk+1)-zb(kk))*AU**2) &
                                    *exp(-a_ext*(x_star-xb(ii))*AU)
            enddo
            image_rt(jj, kk, inu) = I_ray*(yb(jj+1)-yb(jj))*(zb(kk+1)-zb(kk))*(AU/d_obs)**2 !flux per pixel = I * solid angle
        enddo
    enddo
enddo
!$omp end parallel do
end subroutine

subroutine parse_command_line()
!usage: ./mc3d [imagert | imagemc | mctherm | imagetherm] [inu i1 i2 ...] [nphot N]
!   imagert - n_case = 1, source function method + ray-traced images (default)
!   imagemc - n_case = 2, monte carlo picture method
!   mctherm - n_case = 3, dust temperature
!   imagetherm - n_case = 4, dust temperature + imagert (scattered+thermal emission)
!   inu     - wavelength indices to write out (default: all)
!   nphot   - number of photon packets per wavelength, e.g. 100000 or 1e5 (default: value in the code)
implicit none
character(100) :: arg
integer :: iarg, nargs, ival, ios
real(DP) :: rval
logical :: reading_inu, reading_nphot

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
    !images: all wavelengths; thermal modes: only temperatures (every selected wavelength costs nphot packets)
    out_lam = (n_case.eq.1).or.(n_case.eq.2)
    if (n_case.eq.4) call usage('imagetherm needs the image wavelengths, e.g. imagetherm inu 120 250')
else if (.not.any(out_lam)) then
    call usage('no wavelength indices given after "inu"')
endif
select case(n_case)
case(1); print*, 'run mode: imagert (source function method + ray tracing)'
case(2); print*, 'run mode: imagemc (monte carlo picture method)'
case(3); print*, 'run mode: mctherm (dust temperature)'
case(4); print*, 'run mode: imagetherm (dust temperature + imagert (scattered+thermal emission))'
end select
if ((n_case.eq.3).or.(n_case.eq.4)) then
    print*, 'heating packets (wavelengths from the stellar spectrum): ', n_phot, ', per selected wavelength: ', n_phot
else
    print*, 'photon packets per wavelength: ', n_phot
endif
end subroutine

subroutine usage(msg)
implicit none
character(*), intent(in) :: msg
if (len_trim(msg).gt.0) print*, trim(msg)
print*, 'usage: ./mc3d [imagert | imagemc | mctherm | imagetherm] [inu i1 i2 ...] [nphot N]'
print*, '   imagert - source function method + ray-traced images (default)'
print*, '   imagemc - monte carlo picture method'
print*, '   mctherm - dust temperature'
print*, '   imagetherm - dust temperature + imagert (scattered+thermal emission)'
print*, '   inu     - wavelength indices to write out (default: all for images, none for mctherm,'
print*, '             required for imagetherm)'
print*, '   nphot   - photon packets per wavelength (and per time step), e.g. 100000 or 1e5;'
print*, '             mctherm/imagetherm: in addition nphot heating packets in total'
print*, 'with lightcurve.inp in the folder the run is time-dependent'
stop
end subroutine

subroutine deposit(ii, jj, kk, inu_d, heating, val)
!path-length estimator: val = weight*energy*path/(4 pi V) (per time); heating packets add kappa_abs*val
!to the absorbed power, packets of a selected wavelength add val to its J
implicit none
integer, intent(in) :: ii, jj, kk, inu_d
logical, intent(in) :: heating
real(DP), intent(in) :: val
integer :: isel
if (heating) then
    !$omp atomic
    heat_abs(ii,jj,kk) = heat_abs(ii,jj,kk) + kappa_abs(inu_d)*val
else
    isel = isel_of(inu_d)
    if (isel.eq.0) return
    !$omp atomic
    source_func(ii,jj,kk,isel) = source_func(ii,jj,kk,isel) + val
endif
end subroutine

subroutine setup_heating()
!frequency weights for the integrals over the wavelength grid (trapezoidal)
implicit none
integer :: il
if (maxval(kappa_abs).le.0) print*, 'NOTE: kappa_abs = 0 at all wavelengths, the dust is not heated'
if (nlam.lt.2) stop 'mctherm/imagetherm: at least 2 wavelengths are needed in dustkappa.inp'
if (.not.(all(nu(2:nlam).gt.nu(1:nlam-1)).or.all(nu(2:nlam).lt.nu(1:nlam-1)))) &
    stop 'mctherm/imagetherm: the wavelengths in dustkappa.inp must be monotonic'
allocate(dnu_w(nlam), cdf_heat(nlam), p_heat(nlam))
dnu_w(1)    = 0.5*abs(nu(2)-nu(1))
dnu_w(nlam) = 0.5*abs(nu(nlam)-nu(nlam-1))
do il=2, nlam-1
    dnu_w(il) = 0.5*abs(nu(il+1)-nu(il-1))
enddo
end subroutine

subroutine set_heat_cdf(lnu)
!probabilities of the wavelengths for the heating packets: p_i ~ L_nu,i dnu_i
implicit none
real(DP), intent(in) :: lnu(:)
integer :: il
p_heat = lnu*dnu_w
L_heat = sum(p_heat)
if (L_heat.le.0) return
p_heat = p_heat/L_heat
cdf_heat(1) = p_heat(1)
do il=2, nlam
    cdf_heat(il) = cdf_heat(il-1) + p_heat(il)
enddo
cdf_heat(nlam) = 1
end subroutine

function set_heat_cdf_t(t)
!set_heat_cdf for the stellar spectrum at time t of the light curve; returns L_heat
implicit none
real(DP), intent(in) :: t
real(DP) :: set_heat_cdf_t, lnu(nlam)
integer :: il
do il=1, nlam
    lnu(il) = lc_lum(il, t)
enddo
call set_heat_cdf(lnu)
set_heat_cdf_t = L_heat
end function

function draw_heat_lambda()
!wavelength index of a heating packet: smallest i with cdf_heat(i) > xi
implicit none
integer :: draw_heat_lambda, lo, hi, mid
real(DP) :: xi
call random_number(xi)
lo = 0; hi = nlam
do while (hi-lo.gt.1)
    mid = (lo+hi)/2
    if (cdf_heat(mid).gt.xi) then
        hi = mid
    else
        lo = mid
    endif
enddo
draw_heat_lambda = hi
end function

function bytes_packet()
implicit none
real(DP) :: bytes_packet
type(packet) :: p
bytes_packet = storage_size(p)/8
end function

subroutine progress(label, i_done, n_all)
!prints 10, 20, ... 100 % of a loop (called by every iteration, also inside parallel loops)
implicit none
character(*), intent(in) :: label
integer(8), intent(in) :: i_done, n_all
integer(8) :: step
if (n_all.lt.100000) return
step = n_all/10
if (mod(i_done, step).eq.0) print '(3x, a, a, i4, a)', label, ':', int(100*i_done/n_all), ' %'
end subroutine

subroutine memory_check()
!memory of the large arrays, compared with the available memory (Linux: /proc/meminfo)
implicit none
real(DP) :: cells
character(200) :: line
integer :: ios
real(DP) :: kb
cells = real(nx,DP)*ny*nz
mem_static = cells*8*(1 + n_sel)                                   !rho, J of the selected wavelengths
if (thermal) mem_static = mem_static + cells*8*2                   !absorbed power, temperature
if (timedep.and.images_rt) mem_static = mem_static + cells*8*3*n_sel !retarded ray tracing weights and S
mem_avail = 0
open(14, file='/proc/meminfo', status='old', action='read', iostat=ios)
if (ios.eq.0) then
    do
        read(14, '(a)', iostat=ios) line
        if (ios.ne.0) exit
        if (line(1:13).eq.'MemAvailable:') then
            read(line(14:), *, iostat=ios) kb
            if (ios.eq.0) mem_avail = kb*1024
            exit
        endif
    enddo
    close(14)
endif
if (mem_avail.gt.0) then
    print '(a, f8.2, a, f8.2, a)', ' memory: arrays ', mem_static/1e9, ' GB, available ', mem_avail/1e9, ' GB'
    if (mem_static.gt.0.9*mem_avail) then
        print*, 'ERROR: the arrays need more memory than available. Use fewer selected wavelengths (inu)'
        print*, '   or a coarser grid.'
        error stop
    endif
else
    print '(a, f8.2, a)', ' memory: arrays ', mem_static/1e9, ' GB (available memory unknown on this system)'
endif
if (timedep) print '(a, i0, a)', ' (plus ', int(bytes_packet()), ' bytes per packet in flight, shown at every step)'
end subroutine

subroutine setup_grid_density()
!grid.inp + dust_density.inp: read both; neither: default grid with a spherical shell; only one: error
implicit none
logical :: has_grid, has_dens
integer :: ii, jj, kk
real(DP) :: r2

inquire(file=grid_file, exist=has_grid)
inquire(file=dens_file, exist=has_dens)
if (has_grid.and.(.not.has_dens)) call input_help(dens_file//' is missing (found '//grid_file//')')
if (has_dens.and.(.not.has_grid)) call input_help(grid_file//' is missing (found '//dens_file//')')

if (has_grid) then
    call read_grid(grid_file)
else
    !default: (2*n_half_grid+1)**3 cells of 1 AU, centred on the origin
    nx = 2*n_half_grid+1; ny = nx; nz = nx
    allocate(xb(nx+1), yb(ny+1), zb(nz+1))
    do ii=1, nx+1
        xb(ii) = -0.5 + (ii-(n_half_grid+1))
    enddo
    yb = xb; zb = xb
endif
allocate(x(nx), y(ny), z(nz))
x = (xb(1:nx)+xb(2:nx+1))/2 !cell centres [AU]
y = (yb(1:ny)+yb(2:ny+1))/2
z = (zb(1:nz)+zb(2:nz+1))/2
allocate(rho(nx, ny, nz))

if (has_dens) then
    call read_density(dens_file)
    print*, 'grid and dust density read from ', grid_file, ' and ', dens_file, ': ', nx, ' x', ny, ' x', nz, ' cells'
else
    do kk=1, nz
        do jj=1, ny
            do ii=1, nx
                r2 = x(ii)**2+y(jj)**2+z(kk)**2
                if ((r2.lt.r_out**2).and.(r2.gt.r_inn**2)) then
                    rho(ii,jj,kk) = rho0*r_inn**2/r2
                else
                    rho(ii,jj,kk) = 0
                endif
            enddo
        enddo
    enddo
    print*, 'no ', grid_file, '/', dens_file, ': default grid with a spherical shell'
endif

!cell containing the star (on a cell edge: the cell above the edge)
is_star = star_cell(xb, x_star, 'x')
js_star = star_cell(yb, y_star, 'y')
ks_star = star_cell(zb, z_star, 'z')
end subroutine

pure function cell_index(edges, pos, u)
!cell containing pos (binary search); on a cell edge the cell the photon moves into (direction u);
!0 or size(edges): outside the grid / leaving it
implicit none
real(DP), intent(in) :: edges(:), pos, u
integer :: cell_index, lo, hi, mid, ne
ne = size(edges)
if (pos.lt.edges(1)) then
    cell_index = 0
    return
endif
if (pos.gt.edges(ne)) then
    cell_index = ne
    return
endif
lo = 1; hi = ne
do while (hi-lo.gt.1)
    mid = (lo+hi)/2
    if (edges(mid).le.pos) then
        lo = mid
    else
        hi = mid
    endif
enddo
if (edges(hi).le.pos) lo = hi
cell_index = lo
if ((edges(lo).eq.pos).and.(u.le.0)) cell_index = lo - 1
end function

function star_cell(edges, pos, axis)
implicit none
real(DP), intent(in) :: edges(:), pos
character(*), intent(in) :: axis
integer :: star_cell, ii
if ((pos.lt.edges(1)).or.(pos.ge.edges(size(edges)))) then
    print*, 'the star (', axis, ' = ', pos, ' AU) lies outside the grid: ', edges(1), ' .. ', edges(size(edges)), ' AU'
    stop
endif
star_cell = 1
do ii=1, size(edges)-1
    if (edges(ii).le.pos) star_cell = ii
enddo
end function

subroutine read_grid(fname)
!three rows (comment lines starting with # or ! are skipped): cell edges along x, y and z in AU
implicit none
character(*), intent(in) :: fname
character(len=:), allocatable :: line
integer :: iax, n, ios
real(DP), allocatable :: edges(:)
character(1), parameter :: axname(3) = ['x', 'y', 'z']

allocate(character(len=2000000) :: line)
open(12, file=fname, status='old', action='read')
do iax=1, 3
    call read_data_line_or_help(12, line, 'row '//axname(iax)//' of '//fname//' is missing')
    if (len_trim(line).ge.len(line)) call input_help('a row of '//fname//' is too long')
    n = count_values(line)
    if (n.lt.2) call input_help('row '//axname(iax)//' of '//fname//' needs at least 2 cell edges')
    allocate(edges(n))
    read(line, *, iostat=ios) edges
    if (ios.ne.0) call input_help('row '//axname(iax)//' of '//fname//' contains something that is not a number')
    if (any(edges(2:n).le.edges(1:n-1))) call input_help('cell edges in row '//axname(iax)//' of '//fname// &
        ' must increase strictly')
    select case(iax)
    case(1); nx = n-1; call move_alloc(edges, xb)
    case(2); ny = n-1; call move_alloc(edges, yb)
    case(3); nz = n-1; call move_alloc(edges, zb)
    end select
enddo
close(12)
end subroutine

subroutine read_density(fname)
!same layout as dust_temperature3d.out: nz blocks (z from low to high), each with ny lines
!(y from high to low), each line nx values (x from low to high); dust mass density in g cm**-3
implicit none
character(*), intent(in) :: fname
character(len=:), allocatable :: line
character(200) :: msg
integer :: ii, jj, kk, n, ios, iline, ic

allocate(character(len=40*nx+1000) :: line) !generous for nx numbers per line
open(13, file=fname, status='old', action='read')
iline = 0
do kk=1, nz
    do jj=ny, 1, -1
        iline = iline + 1
        write(msg, '(a, i0, a, i0, a, i0, a)') fname//' has fewer data lines than ny*nz = ', ny*nz, &
            ' (stopped at line ', iline, ' = block ', kk, ')'
        call read_data_line_or_help(13, line, trim(msg))
        if (len_trim(line).ge.len(line)) call input_help('a line of '//fname//' is too long (more than nx values?)')
        n = count_values(line)
        if (n.ne.nx) then
            write(msg, '(a, i0, a, i0, a, i0, a)') 'data line ', iline, ' of '//fname//' has ', n, &
                ' values, expected nx = ', nx, ' (number of x cells in grid.inp)'
            call input_help(trim(msg))
        endif
        read(line, *, iostat=ios) (rho(ii,jj,kk), ii=1, nx)
        if (ios.ne.0) then
            write(msg, '(a, i0, a)') 'data line ', iline, ' of '//fname//' contains something that is not a number'
            call input_help(trim(msg))
        endif
    enddo
enddo
!nothing but comments or blank lines may follow
do
    read(13, '(a)', iostat=ios) line
    if (ios.ne.0) exit
    ic = verify(line, ' '//char(9))
    if (ic.eq.0) cycle
    if ((line(ic:ic).eq.'#').or.(line(ic:ic).eq.'!')) cycle
    write(msg, '(a, i0, a)') fname//' has more data lines than ny*nz = ', ny*nz, ' (grid.inp and the density do not match)'
    call input_help(trim(msg))
enddo
close(13)
if (any(rho.ne.rho)) call input_help(fname//' contains NaN values')
if (any(rho.lt.0)) call input_help(fname//' contains negative densities')
end subroutine

subroutine read_data_line_or_help(unit, line, msg)
!next line that is not blank and not a comment; at the end of the file: input_help(msg)
implicit none
integer, intent(in) :: unit
character(*), intent(out) :: line
character(*), intent(in) :: msg
integer :: ios, ic
do
    read(unit, '(a)', iostat=ios) line
    if (ios.ne.0) call input_help(msg)
    ic = verify(line, ' '//char(9)) !first character that is not blank
    if (ic.eq.0) cycle
    if ((line(ic:ic).eq.'#').or.(line(ic:ic).eq.'!')) cycle
    exit
enddo
end subroutine

function count_values(line)
!number of entries separated by blanks, tabs or commas
implicit none
character(*), intent(in) :: line
integer :: count_values, ic
logical :: in_value, sep
count_values = 0
in_value = .false.
do ic=1, len_trim(line)
    sep = (line(ic:ic).eq.' ').or.(line(ic:ic).eq.',').or.(line(ic:ic).eq.char(9))
    if ((.not.sep).and.(.not.in_value)) count_values = count_values + 1
    in_value = .not.sep
enddo
end function

subroutine input_help(msg)
!error in the model input: explain both files and stop
implicit none
character(*), intent(in) :: msg
print*, 'ERROR: ', trim(msg)
print*, ''
print*, 'A model on your own grid needs both files in the run folder (without them the default spherical shell is used):'
print*, ''
print*, 'grid.inp - cell edges in AU, three rows: x edges, y edges, z edges (strictly increasing,'
print*, '   at least 2 per row; nx = number of x edges - 1 etc.; the observer looks along -x;'
print*, '   the star at (x_star, y_star, z_star) must lie inside the grid). Example, 3 x 2 x 2 cells:'
print*, '      -1.5 -0.5 0.5 1.5'
print*, '      -1 0 1'
print*, '      -1 0 1'
print*, ''
print*, 'dust_density.inp - dust mass density in g cm**-3, same layout as dust_temperature3d.out:'
print*, '   nz blocks (z from low to high), each block ny lines (y from HIGH to low),'
print*, '   each line nx values (x from low to high), i.e. ny*nz lines of nx values.'
print*, '   From a numpy array rho[x, y, z]:'
print*, '      np.savetxt("dust_density.inp", rho[:, ::-1, :].transpose(2, 1, 0).reshape(-1, rho.shape[0]))'
print*, ''
print*, 'Lines starting with # or ! are comments in both files.'
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
!int kappa_abs J_nu dnu is accumulated directly in heat_abs by the heating packets (direct + scattered starlight).
!Heating by the dust's own thermal emission is not included (dust optically thin in the IR).
!The CMB (T_cmb) is added as an isotropic background to the absorbed energy, so T >= T_cmb.
!The frequency integrals use trapezoidal weights on the wavelength grid of the opacity file,
!which must be monotonic and cover both the stellar spectrum and the dust emission.
implicit none
integer, parameter :: n_T = 2000
real(DP), parameter :: T_min = 1, T_max = 1e5 ![T] = K, range of the emission table
real(DP) :: T_tab(n_T), E_tab(n_T), E_abs, E_cmb
integer :: it, ii, jj, kk, lo, hi, mid

temp = 0
if (maxval(kappa_abs).le.0) then
    print*, 'dust_temp: kappa_abs = 0 at all wavelengths, no dust heating, temperature set to 0'
    return
endif
!emission table E(T) = sum kappa_abs B(T) dnu, log-spaced in T
do it=1, n_T
    T_tab(it) = T_min*(T_max/T_min)**(real(it-1,DP)/(n_T-1))
    E_tab(it) = sum(kappa_abs*B(nu, T_tab(it))*dnu_w)
enddo
E_cmb = sum(kappa_abs*B(nu, T_cmb)*dnu_w)

!$omp parallel do collapse(3) private(E_abs, lo, hi, mid) schedule(static)
do kk=1, nz
    do jj=1, ny
        do ii=1, nx
            if (rho(ii,jj,kk).le.0) cycle !no dust -> temperature left at 0
            E_abs = heat_abs(ii,jj,kk) + E_cmb
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

!=====================================================================================
! time-dependent radiative transfer (used when lightcurve.inp is present)
!=====================================================================================

subroutine read_lightcurve(fname)
!lightcurve.inp, lines starting with '#' or '!' are comments:
!   iformat     (1: L_star [Lsun], spectrum shape of the code T_star   2: R_star [Rsun], T_star [K])
!   nt          (number of rows that follow)
!   itime  time[s]  L_star            (iformat = 1)
!   itime  time[s]  R_star  T_star    (iformat = 2)
!itime starts at 0: the run starts from the equilibrium radiation field for the first row
!itime starts at 1: the run starts at time 0 from a dark cloud (dust at T_cmb), first time must be > 0
implicit none
character(*), intent(in) :: fname
character(512) :: line
integer :: ir

open(11, file=fname, status='old', action='read')
call read_data_line(11, line); read(line, *) lc_iformat
if ((lc_iformat.ne.1).and.(lc_iformat.ne.2)) stop 'read_lightcurve: iformat must be 1 or 2'
call read_data_line(11, line); read(line, *) lc_n
if (lc_n.lt.1) stop 'read_lightcurve: number of time steps must be >= 1'
allocate(lc_it(lc_n), lc_t(lc_n), lc_Lb(lc_n), lc_Rs(lc_n), lc_Ts(lc_n))
lc_Lb = 0; lc_Rs = 0; lc_Ts = 0
do ir=1, lc_n
    call read_data_line(11, line)
    if (lc_iformat.eq.1) then
        read(line, *) lc_it(ir), lc_t(ir), lc_Lb(ir)
        if (lc_Lb(ir).lt.0) stop 'read_lightcurve: L_star must be >= 0'
    else
        read(line, *) lc_it(ir), lc_t(ir), lc_Rs(ir), lc_Ts(ir)
        if ((lc_Rs(ir).lt.0).or.(lc_Ts(ir).lt.0)) stop 'read_lightcurve: R_star and T_star must be >= 0'
    endif
enddo
close(11)

if ((lc_it(1).ne.0).and.(lc_it(1).ne.1)) stop 'read_lightcurve: itime must start with 0 or 1'
do ir=2, lc_n
    if (lc_it(ir).ne.lc_it(ir-1)+1) stop 'read_lightcurve: itime must increase by 1 from row to row'
    if (lc_t(ir).le.lc_t(ir-1)) stop 'read_lightcurve: time must increase strictly'
enddo
eq_start = (lc_it(1).eq.0)
if (.not.eq_start) then
    if (lc_t(1).eq.0) stop 'read_lightcurve: itime = 1 cannot have time = 0 (the dark start is at time 0)'
    if (lc_t(1).lt.0) stop 'read_lightcurve: itime = 1 must have time > 0 (the dark start is at time 0)'
endif
if (eq_start) then
    r_first = 2  !row 1 is the equilibrium state, steps start from row 2
    t_init = lc_t(1)
else
    r_first = 1
    t_init = 0
endif
end subroutine

function lc_lum(inu, t)
!stellar L_nu [erg s-1 Hz-1] at time t, linear interpolation between the rows of lightcurve.inp,
!constant before the first and after the last row
implicit none
integer, intent(in) :: inu
real(DP), intent(in) :: t
real(DP) :: lc_lum, f, Lb, Rs, Ts
integer :: ir, r1, r2

if (t.le.lc_t(1)) then
    r1 = 1; r2 = 1; f = 0
else if (t.ge.lc_t(lc_n)) then
    r1 = lc_n; r2 = lc_n; f = 0
else
    r2 = lc_n
    do ir=2, lc_n
        if (t.le.lc_t(ir)) then
            r2 = ir
            exit
        endif
    enddo
    r1 = r2 - 1
    f = (t - lc_t(r1))/(lc_t(r2) - lc_t(r1))
endif
if (lc_iformat.eq.1) then
    Lb = (1-f)*lc_Lb(r1) + f*lc_Lb(r2)
    lc_lum = L_star(inu)*Lb*LSun/(4*pi*R_star**2*sigma_SB*T_star**4) !shape of B(T_star), scaled to Lb
else
    Rs = ((1-f)*lc_Rs(r1) + f*lc_Rs(r2))*RSun
    Ts = (1-f)*lc_Ts(r1) + f*lc_Ts(r2)
    lc_lum = 4*pi*Rs**2*pi*B(nu(inu), Ts)
endif
end function

subroutine run_timedep()
!time stepping with a census of packets in flight (option B):
!   every step: propagate the packets in flight to the end of the step, emit new packets during the step,
!   keep the survivors for the next step, J is the average over the step
!   images: imagert/imagetherm - S of every step is added to the images in retarded time
!                                (step averages assigned to the step middles, linear interpolation in between),
!           imagemc - escaping packets are binned in observer time
!   after the last row the run continues (light curve kept flat) until all light that
!   reaches the observer before the last time moment has left the grid
implicit none
integer :: ir, isel, m, id, n_d, inu_p, f_eq
real(DP) :: pix_area
integer(8) :: n_snap, n_snap_old, n_eq
real(DP) :: t_prev, t_last, t_end, dt_d, age, xi
type(packet) :: p
logical :: snapped

!packets per time step: nphot per selected wavelength, plus nphot heating packets in the thermal modes
n_act = n_sel
allocate(act(n_act)); act = sel
n_tot = int(n_phot,8)*n_act

print*, 'lightcurve.inp found: time-dependent run with ', lc_n, ' time moments'
if (eq_start) then
    print*, 'start: equilibrium radiation field at itime = 0, time = ', lc_t(1), ' s'
else
    print*, 'start: dark cloud at time = 0, dust temperature ', T_cmb, ' K'
endif
print*, 'packets per time step: ', n_tot, ' for the selected wavelengths', merge(int(n_phot,8), 0_8, thermal), ' for heating'

if (n_sel.gt.0) call execute_command_line('mkdir -p J_timeseries')
if (thermal) call execute_command_line('mkdir -p T_timeseries')
if (images_rt.or.images_mc) then
    allocate(img_t(ny, nz, lc_n, n_sel))
    img_t = 0
endif
if (images_rt) then
    call retarded_weights()
    allocate(s_prev(nx, ny, nz, n_sel), s_cur(nx, ny, nz, n_sel))
    !dark start: no scattered light before time 0
    s_prev = 0
    c_prev = t_init
endif
if (images_mc) then
    !observer time bins: row r collects (tbin_lo(r), tbin_hi(r)], the equilibrium row is the steady image
    allocate(tbin_lo(lc_n), tbin_hi(lc_n))
    do ir=1, lc_n
        tbin_hi(ir) = lc_t(ir)
        if (ir.eq.1) then
            tbin_lo(ir) = t_init
        else
            tbin_lo(ir) = lc_t(ir-1)
        endif
    enddo
endif
allocate(pk_new(1024))
n_new = 0

!----- equilibrium start: steady-state packets with random ages, snapshot at t_init gives the census
if (eq_start) then
    t_burn = n_burn*sqrt((xb(nx+1)-xb(1))**2+(yb(ny+1)-yb(1))**2+(zb(nz+1)-zb(1))**2)*AU/clight !grid diagonal
    !as many packets per unit time as in the first step, so the census at itime = 0 is as good as later ones
    if (lc_n.ge.2) then
        f_eq = int(min(max(1.0_DP, t_burn/(lc_t(2)-lc_t(1))), real(f_eq_max,DP)))
    else
        f_eq = 1
    endif
    n_eq = n_tot*f_eq
    print*, 'equilibrium run with ', f_eq, ' x the packets of a time step'
    if (f_eq.eq.f_eq_max) print*, 'NOTE: first time step is short, packets in flight at itime = 0 are noisier than later'
    inv_dt = 1/t_burn
    source_func = 0
    if (thermal) heat_abs = 0
    n_snap = 0; n_snap_old = 0
    n_long_walk = 0
    !packets of the selected wavelengths
    !$omp parallel do schedule(dynamic, 100) private(p, snapped, age, xi, inu_p) reduction(+:n_snap, n_snap_old)
    do ip=1, n_eq
        inu_p = act(int(mod(ip-1, int(n_act,8)))+1)
        call random_number(xi)
        age = xi*t_burn
        !steady power L/(n_phot*f_eq) carried as energy over t_burn
        call emit_packet(p, inu_p, t_init-age, lc_lum(inu_p, t_init)*t_burn/(n_phot*f_eq), .false.)
        call propagate(p, t_init, .true., snapped)
        if (snapped) then
            n_snap = n_snap + 1
            if (age.gt.t_burn/2) n_snap_old = n_snap_old + 1
        endif
        call progress('equilibrium, selected wavelengths', ip, n_eq)
    enddo
    !$omp end parallel do
    !heating packets: wavelengths from the stellar spectrum at itime = 0
    !(nested if: Fortran does not guarantee short-circuit evaluation of .and.)
    L_heat = 0
    if (thermal) L_heat = set_heat_cdf_t(t_init)
    if (L_heat.gt.0) then
        !$omp parallel do schedule(dynamic, 100) private(p, snapped, age, xi, inu_p) reduction(+:n_snap, n_snap_old)
        do ip=1, int(n_phot,8)*f_eq
            inu_p = draw_heat_lambda()
            call random_number(xi)
            age = xi*t_burn
            call emit_packet(p, inu_p, t_init-age, &
                lc_lum(inu_p, t_init)*dnu_w(inu_p)*t_burn/(int(n_phot,8)*f_eq*p_heat(inu_p)), .true.)
            call propagate(p, t_init, .true., snapped)
            if (snapped) then
                n_snap = n_snap + 1
                if (age.gt.t_burn/2) n_snap_old = n_snap_old + 1
            endif
            call progress('equilibrium, heating packets', ip, int(n_phot,8)*f_eq)
        enddo
        !$omp end parallel do
    endif
    print*, 'equilibrium done, packets in flight at itime = 0: ', n_new, ' (', n_new*bytes_packet()/1e6, ' MB)'
    if (n_snap.gt.0) then
        if (real(n_snap_old,DP)/n_snap.gt.0.01) print*, 'WARNING: ', 100*real(n_snap_old,DP)/n_snap, &
            ' % of packets in flight are older than t_burn/2, increase n_burn'
    endif
    if (n_long_walk.gt.0) then
        print*, 'NOTE: ', n_long_walk, ' random walks in optically thick cells last longer than t_burn = ', t_burn, ' s.'
        print*, '   The radiation trapped in these cells is in the equilibrium J of itime 0 but not in the packets'
        print*, '   in flight, so J inside them drops after itime 0 and builds up again on the diffusion time.'
        print*, '   Images and optically thin regions are not affected; with absorption this is rare.'
    endif
    call post_step(lc_it(1), -huge(1.0_DP), t_init, .true.)
endif

!----- time steps between the rows of lightcurve.inp
t_prev = t_init
do ir=r_first, lc_n
    call time_step(t_prev, lc_t(ir), .true.)
    print '(a, i6, a, es12.4, a, i12, a, f10.1, a)', ' itime ', lc_it(ir), '  time [s] ', lc_t(ir), &
        '  packets in flight ', n_new, ' (', n_new*bytes_packet()/1e6, ' MB)'
    call post_step(lc_it(ir), t_prev, lc_t(ir), .true.)
    t_prev = lc_t(ir)
enddo

!----- drain: light emitted before the last time moment still has to reach the images; the light curve is
!continued flat (light emitted later arrives at the observer after the last time moment, but a flat
!continuation keeps the step averages of J around the last moment unbiased)
if (images_rt.or.images_mc) then
    t_last = lc_t(lc_n)
    t_end = t_last + (xb(nx+1)-xb(1))*AU/clight
    if (lc_n.ge.2) then
        dt_d = max(lc_t(lc_n)-lc_t(lc_n-1), minval(xb(2:nx+1)-xb(1:nx))*AU/clight)
    else if (eq_start) then
        dt_d = t_end - t_last !only the equilibrium: one step is enough
    else
        dt_d = max(lc_t(lc_n)-t_init, minval(xb(2:nx+1)-xb(1:nx))*AU/clight)
    endif
    n_d = max(1, ceiling((t_end-t_last)/dt_d))
    dt_d = (t_end-t_last)/n_d
    print*, 'continuing for ', n_d, ' more steps until the light reaching the last image has left the grid'
    do id=1, n_d
        call time_step(t_prev, t_prev+dt_d, .true.)
        call post_step(-1, t_prev, t_prev+dt_d, .false.)
        t_prev = t_prev + dt_d
    enddo
    !retarded times after the middle of the last step: source function of the last step
    if (images_rt) call retarded_add(c_prev, huge(1.0_DP), c_prev, c_prev, .false.)
endif

!----- images, flux per pixel [erg s-1 cm-2 Hz-1], rows k = N..1, columns j = 1..N
if (images_rt.or.images_mc) then
    do isel=1, n_sel
        inu = sel(isel)
        do m=1, lc_n
            if (images_rt) then
                !unresolved star at the observer time of the image
                pix_area = (yb(js_star+1)-yb(js_star))*(zb(ks_star+1)-zb(ks_star))*AU**2
                img_t(js_star, ks_star, m, isel) = img_t(js_star, ks_star, m, isel) &
                    + lc_lum(inu, lc_t(m))/(4*pi*pix_area)*star_att(isel)
                !flux per pixel = intensity * solid angle of the pixel
                do k=1, nz
                    do j=1, ny
                        img_t(j,k,m,isel) = img_t(j,k,m,isel)*(yb(j+1)-yb(j))*(zb(k+1)-zb(k))*(AU/d_obs)**2
                    enddo
                enddo
            endif
            select case(n_case)
            case(1); write(filename, '(a, i4.4, a, i0, a)') 'image_raytrace_', inu, '_', lc_it(m), '.out'
            case(2); write(filename, '(a, i4.4, a, i0, a)') 'image_montecarlo_', inu, '_', lc_it(m), '.out'
            case(4); write(filename, '(a, i4.4, a, i0, a)') 'image_raytrace_therm_', inu, '_', lc_it(m), '.out'
            end select
            open(430, file=filename)
            do k=nz, 1, -1
                write(430, '(*(es13.3e3))') (img_t(j, k, m, isel), j=1, ny)
            enddo
            close(430)
        enddo
    enddo
endif
end subroutine

subroutine time_step(ta, tb, emit)
!propagate the census from ta to tb, emit new packets uniformly in (ta, tb] if emit
implicit none
real(DP), intent(in) :: ta, tb
logical, intent(in) :: emit
type(packet) :: p
integer(8) :: ipk, n_old
real(DP) :: xi, te
integer :: inu_p
logical :: snapped

inv_dt = 1/(tb-ta)
source_func = 0
if (thermal) heat_abs = 0
call move_alloc(pk_new, pk_old)
n_old = n_new
n_old_census = n_old
allocate(pk_new(max(1024_8, n_old)))
n_new = 0
!$omp parallel do schedule(dynamic, 100) private(p, snapped)
do ipk=1, n_old
    p = pk_old(ipk)
    call propagate(p, tb, .false., snapped)
enddo
!$omp end parallel do
deallocate(pk_old)
n_old_census = 0
if (emit) then
    !$omp parallel do schedule(dynamic, 100) private(p, snapped, xi, te, inu_p)
    do ipk=1, n_tot
        inu_p = act(int(mod(ipk-1, int(n_act,8)))+1)
        call random_number(xi)
        te = ta + xi*(tb-ta)
        !energy per packet: expectation over te gives int L dt / n_phot
        call emit_packet(p, inu_p, te, lc_lum(inu_p, te)*(tb-ta)/n_phot, .false.)
        call propagate(p, tb, .false., snapped)
        if (n_tot.ge.2000000) call progress('time step, selected wavelengths', ipk, n_tot)
    enddo
    !$omp end parallel do
    !heating packets: wavelength drawn from the spectrum in the middle of the step, energy weighted
    !with L_nu(te)/p_nu so that the expectation is int L_nu dnu dt / n_phot for any te
    L_heat = 0
    if (thermal) L_heat = set_heat_cdf_t((ta+tb)/2)
    if (L_heat.gt.0) then
        !$omp parallel do schedule(dynamic, 100) private(p, snapped, xi, te, inu_p)
        do ipk=1, int(n_phot,8)
            inu_p = draw_heat_lambda()
            call random_number(xi)
            te = ta + xi*(tb-ta)
            call emit_packet(p, inu_p, te, lc_lum(inu_p, te)*dnu_w(inu_p)*(tb-ta)/(n_phot*p_heat(inu_p)), .true.)
            call propagate(p, tb, .false., snapped)
            if (n_phot.ge.2000000) call progress('time step, heating packets', ipk, int(n_phot,8))
        enddo
        !$omp end parallel do
    endif
endif
end subroutine

subroutine post_step(itime, ta, tb, dump)
!source_func holds J averaged over (ta, tb]
implicit none
integer, intent(in) :: itime
real(DP), intent(in) :: ta, tb
logical, intent(in) :: dump
character(100) :: fname
integer :: ii, jj, kk, il

if (dump.and.(n_sel.gt.0)) then
    !J of the selected wavelengths: for each wavelength the full cube, blocks of constant k, rows j = N..1, columns i = 1..N
    write(fname, '(a, i0, a)') 'J_timeseries/J_', itime, '.out'
    open(420, file=fname)
    do il=1, n_sel
        do kk=1, nz
            do jj=ny, 1, -1
                write(420, '(*(es13.3e3))') (source_func(ii,jj,kk,il), ii=1, nx)
            enddo
        enddo
    enddo
    close(420)
endif
if (thermal) then
    call dust_temp()
    if (dump) then
        write(fname, '(a, i0, a)') 'T_timeseries/T_', itime, '.out'
        open(421, file=fname)
        do kk=1, nz
            do jj=ny, 1, -1
                write(421, '(*(es13.3e3))') (temp(ii,jj,kk), ii=1, nx)
            enddo
        enddo
        close(421)
    endif
endif
if (images_rt) then
    call source_selected()
    if (ta.eq.-huge(1.0_DP)) then
        !equilibrium: constant source function for all retarded times up to tb
        s_prev = s_cur
        call retarded_add(-huge(1.0_DP), tb, tb, tb, .false.)
        c_prev = tb
    else
        !step average assigned to the step middle, linear interpolation from the previous step
        call retarded_add(c_prev, (ta+tb)/2, c_prev, (ta+tb)/2, .true.)
        s_prev = s_cur
        c_prev = (ta+tb)/2
    endif
endif
end subroutine

subroutine source_selected()
!source function of the selected wavelengths from J (and T) of the current step
implicit none
integer :: isel, inu_s, ii, jj, kk
!$omp parallel do collapse(2) private(inu_s, ii) schedule(static)
do isel=1, n_sel
    do kk=1, nz
        inu_s = sel(isel)
        do jj=1, ny
            do ii=1, nx
                s_cur(ii,jj,kk,isel) = albedo(inu_s)*source_func(ii,jj,kk,isel)
                if (n_case.eq.4) s_cur(ii,jj,kk,isel) = s_cur(ii,jj,kk,isel) &
                    + (1-albedo(inu_s))*B(nu(inu_s), temp(ii,jj,kk))
            enddo
        enddo
    enddo
enddo
!$omp end parallel do
end subroutine

subroutine emit_packet(p, inu_p, te, energy, heating)
!new packet at the star with energy [erg Hz-1], or [erg] for a heating packet
implicit none
type(packet), intent(out) :: p
integer, intent(in) :: inu_p
real(DP), intent(in) :: te, energy
logical, intent(in) :: heating
p%x = x_star; p%y = y_star; p%z = z_star
p%ys = y_star; p%zs = z_star
p%t = te; p%tobs = te
p%w = 1
p%l_walk = 0; p%r_walk = 0
p%heat = heating
p%e = energy
p%inu = inu_p
call new_direction(p)
end subroutine

subroutine new_direction(p)
!isotropic direction and optical depth to the next interaction
implicit none
type(packet), intent(inout) :: p
real(DP) :: xi, xi1, xi2, phi, cos_theta, sin_theta
call random_number(xi)
p%tau = -log(xi)
call random_number(xi1)
call random_number(xi2)
phi = 2*pi*xi1
cos_theta = 2*xi2 - 1
sin_theta = sqrt(1 - cos_theta**2)
p%ux = sin_theta*cos(phi)
p%uy = sin_theta*sin(phi)
p%uz = cos_theta
!current cell for the new flight (on a wall: the cell the photon moves into)
p%ix = cell_index(xb, p%x, p%ux)
p%iy = cell_index(yb, p%y, p%uy)
p%iz = cell_index(zb, p%z, p%uz)
end subroutine

subroutine append_census(p)
!add a packet in flight to the census of the next step (grows the array if needed)
implicit none
type(packet), intent(in) :: p
type(packet), allocatable :: tmp(:)
!$omp critical (census_append)
if (n_new.ge.size(pk_new, kind=8)) then
    !growing needs the old and the new array at the same time (plus the census of the previous step)
    if (mem_avail.gt.0) then
        if (mem_static + (3*size(pk_new, kind=8) + n_old_census)*bytes_packet().gt.0.9*mem_avail) then
            print*, 'ERROR: the packets in flight would need more memory than available: ', &
                (3*size(pk_new, kind=8) + n_old_census)*bytes_packet()/1e9, ' GB for ', n_new, ' packets'
            print*, '   reduce nphot, use a longer first time step (equilibrium start), or fewer selected wavelengths'
            error stop
        endif
    endif
    allocate(tmp(max(2*size(pk_new, kind=8), 1024_8)))
    tmp(1:n_new) = pk_new(1:n_new)
    call move_alloc(tmp, pk_new)
endif
n_new = n_new + 1
pk_new(n_new) = p
!$omp end critical (census_append)
end subroutine

subroutine propagate(p, t_stop, steady, snapped)
!move packet p until it escapes, is killed by russian roulette or reaches t_stop
!steady = .true. (equilibrium run): at t_stop a copy goes to the census and the packet continues until it escapes
!steady = .false. (time step): at t_stop the packet goes to the census and propagate returns
implicit none
type(packet), intent(inout) :: p
real(DP), intent(in) :: t_stop
logical, intent(in) :: steady
logical, intent(out) :: snapped
real(DP) :: xw, yw, zw, sx, sy, sz, s_wall, s_int, s_t, s, a_ext, vol_c, xi, x1, y1, z1
real(DP) :: r0, seg, a_abs, ex, ey, ez !modified random walk
integer :: im(1), jm(1), km(1), event
logical :: stop_check

snapped = .false.
stop_check = .true.
do
    !current cell, kept in the packet (search in new_direction, updated at every wall crossing)
    if ((p%ix.lt.1).or.(p%ix.gt.nx).or.(p%iy.lt.1).or.(p%iy.gt.ny).or.(p%iz.lt.1).or.(p%iz.gt.nz)) then
        if ((p%ix.lt.1).and.images_mc.and.(p%ux+1.le.dmu_obs)) call tally_mc(p, steady)
        return
    endif
    im(1) = p%ix; jm(1) = p%iy; km(1) = p%iz

    !modified random walk in optically thick cells: start a walk, or continue one from the previous time step
    a_ext = rho(im(1), jm(1), km(1))*kappa_ext(p%inu)
    if (p%l_walk.le.0) then
        r0 = min(p%x-xb(im(1)), xb(im(1)+1)-p%x, p%y-yb(jm(1)), yb(jm(1)+1)-p%y, p%z-zb(km(1)), zb(km(1)+1)-p%z)
        if (a_ext*r0*AU.gt.mrw_gamma) then
            p%l_walk = mrw_length(a_ext, r0*AU)
            p%r_walk = r0
            if (steady.and.(p%l_walk/clight.gt.t_burn)) then
                !$omp atomic
                n_long_walk = n_long_walk + 1
            endif
        endif
    endif
    if (p%l_walk.gt.0) then
        !the walk takes l_walk/c; it may end at the end of the time step
        seg = p%l_walk; event = 0
        if (stop_check) then
            s_t = max((t_stop-p%t)*clight, 0.0_DP)
            if (s_t.lt.seg) then
                seg = s_t; event = 2
            endif
        endif
        a_abs = rho(im(1), jm(1), km(1))*kappa_abs(p%inu)
        vol_c = abs(xb(im(1)+1)-xb(im(1)))*abs(yb(jm(1)+1)-yb(jm(1)))*abs(zb(km(1)+1)-zb(km(1)))*AU**3
        call deposit(im(1), jm(1), km(1), p%inu, p%heat, p%e*p%w*mrw_eff_length(a_abs, seg)*inv_dt/4./pi/vol_c)
        p%w = p%w*exp(-a_abs*seg)
        p%t = p%t + seg/clight
        if (event.eq.2) then
            !still inside the sphere at the end of the step: the rest of the walk is done in the next step
            p%l_walk = p%l_walk - seg
            p%t = t_stop
            snapped = .true.
            call append_census(p)
            if (.not.steady) return
            stop_check = .false.
            cycle
        endif
        !walk finished: the photon leaves the sphere surface outwards (cosine-weighted around the normal)
        p%l_walk = 0
        call random_dir(ex, ey, ez)
        p%x = p%x + p%r_walk*ex; p%y = p%y + p%r_walk*ey; p%z = p%z + p%r_walk*ez
        p%x = min(max(p%x, xb(im(1))), xb(im(1)+1)) !the sphere lies inside the cell
        p%y = min(max(p%y, yb(jm(1))), yb(jm(1)+1))
        p%z = min(max(p%z, zb(km(1))), zb(km(1)+1))
        if (p%w.lt.w_min) then
            call random_number(xi)
            if (xi.ge.p_survive) return
            p%w = p%w/p_survive
        endif
        p%ys = p%y; p%zs = p%z
        p%tobs = p%t + (p%x-x_star)*AU/clight
        call outward_dir(ex, ey, ez, p%ux, p%uy, p%uz)
        call random_number(xi)
        p%tau = -log(xi)
        cycle
    endif

    !path length (in AU) to the x-, y- and z-walls of the current cell
    if (p%ux.gt.0) then
        xw = xb(im(1)+1); sx = max((xw-p%x)/p%ux, 0.0_DP)
    else if (p%ux.lt.0) then
        xw = xb(im(1));   sx = max((xw-p%x)/p%ux, 0.0_DP)
    else
        xw = p%x;         sx = huge(1.0_DP)
    endif
    if (p%uy.gt.0) then
        yw = yb(jm(1)+1); sy = max((yw-p%y)/p%uy, 0.0_DP)
    else if (p%uy.lt.0) then
        yw = yb(jm(1));   sy = max((yw-p%y)/p%uy, 0.0_DP)
    else
        yw = p%y;         sy = huge(1.0_DP)
    endif
    if (p%uz.gt.0) then
        zw = zb(km(1)+1); sz = max((zw-p%z)/p%uz, 0.0_DP)
    else if (p%uz.lt.0) then
        zw = zb(km(1));   sz = max((zw-p%z)/p%uz, 0.0_DP)
    else
        zw = p%z;         sz = huge(1.0_DP)
    endif
    s_wall = min(sx, sy, sz)

    !closest event: wall (0), interaction (1) or end of the time step (2)
    if (a_ext.gt.0) then
        s_int = p%tau/(a_ext*AU)
    else
        s_int = huge(1.0_DP)
    endif
    s = s_wall; event = 0
    if (s_int.lt.s) then
        s = s_int; event = 1
    endif
    if (stop_check) then
        s_t = max((t_stop-p%t)*clight/AU, 0.0_DP)
        if (s_t.le.s) then
            s = s_t; event = 2
        endif
    endif

    !path-length estimator of J averaged over the time step (or over t_burn in the equilibrium run)
    vol_c = abs(xb(im(1)+1)-xb(im(1)))*abs(yb(jm(1)+1)-yb(jm(1)))*abs(zb(km(1)+1)-zb(km(1)))*AU**3
    call deposit(im(1), jm(1), km(1), p%inu, p%heat, p%e*p%w*s*AU*inv_dt/4./pi/vol_c)

    if (event.eq.0) then
        x1 = p%x + s*p%ux
        y1 = p%y + s*p%uy
        z1 = p%z + s*p%uz
        if (sx.eq.s) x1 = xw !put the photon exactly on the crossed wall
        if (sy.eq.s) y1 = yw
        if (sz.eq.s) z1 = zw
        p%x = x1; p%y = y1; p%z = z1
        !wall crossing: neighbouring cell
        if (sx.eq.s) p%ix = p%ix + merge(1, -1, p%ux.gt.0)
        if (sy.eq.s) p%iy = p%iy + merge(1, -1, p%uy.gt.0)
        if (sz.eq.s) p%iz = p%iz + merge(1, -1, p%uz.gt.0)
    else
        p%x = p%x + s*p%ux
        p%y = p%y + s*p%uy
        p%z = p%z + s*p%uz
    endif
    p%t = p%t + s*AU/clight
    p%tau = p%tau - a_ext*s*AU

    if (event.eq.2) then
        p%t = t_stop
        snapped = .true.
        call append_census(p)
        if (.not.steady) return
        stop_check = .false.
    else if (event.eq.1) then
        !interaction: absorbed fraction removed from the packet, the rest is scattered
        p%w = p%w*albedo(p%inu)
        if (p%w.lt.w_min) then
            call random_number(xi)
            if (xi.ge.p_survive) return
            p%w = p%w/p_survive
        endif
        !observer time of light leaving this point towards -x, and its position on the sky
        p%ys = p%y; p%zs = p%z
        p%tobs = p%t + (p%x-x_star)*AU/clight
        call new_direction(p)
    endif
enddo
end subroutine

subroutine tally_mc(p, steady)
!packet leaving through the -x face inside the cone towards the observer
implicit none
type(packet), intent(in) :: p
logical, intent(in) :: steady
integer :: jm(1), km(1), m, ir, isel
real(DP) :: fac

if (p%heat) return
isel = isel_of(p%inu)
if (isel.eq.0) return
jm(1) = cell_index(yb, p%ys, -1.0_DP)
km(1) = cell_index(zb, p%zs, -1.0_DP)
if ((jm(1).lt.1).or.(jm(1).gt.ny).or.(km(1).lt.1).or.(km(1).gt.nz)) return
if (steady) then
    m = 1 !equilibrium image
    fac = inv_dt
else
    m = 0
    do ir=r_first, lc_n
        if ((p%tobs.gt.tbin_lo(ir)).and.(p%tobs.le.tbin_hi(ir))) then
            m = ir
            exit
        endif
    enddo
    if (m.eq.0) return
    fac = 1/(tbin_hi(m)-tbin_lo(m)) !flux averaged over the observer time bin
endif
!energy per steradian of the cone (solid angle 2*pi*dmu_obs) per time, divided by d_obs**2
!$omp atomic
img_t(jm(1), km(1), m, isel) = img_t(jm(1), km(1), m, isel) + p%e*p%w*fac/(2*pi*dmu_obs)/d_obs**2
end subroutine

subroutine retarded_weights()
!static weights of the formal solution along -x: I = sum_i S_i (1-exp(-dtau_i)) exp(-tau in front of cell i)
implicit none
integer :: isel, ii, jj, kk, inu_s
real(DP) :: dtau, tau_f

allocate(w_rt(nx, ny, nz, n_sel), star_att(n_sel), delta_x(nx))
do ii=1, nx
    delta_x(ii) = (x(ii)-x_star)*AU/clight !extra light travel time to the observer compared with the star
enddo
do isel=1, n_sel
    inu_s = sel(isel)
    do kk=1, nz
        do jj=1, ny
            tau_f = 0
            do ii=1, nx
                dtau = rho(ii,jj,kk)*kappa_ext(inu_s)*(xb(ii+1)-xb(ii))*AU
                w_rt(ii,jj,kk,isel) = (1-exp(-dtau))*exp(-tau_f)
                if ((kk.eq.ks_star).and.(jj.eq.js_star).and.(ii.eq.is_star)) &
                    star_att(isel) = exp(-rho(ii,jj,kk)*kappa_ext(inu_s)*(x_star-xb(ii))*AU)*exp(-tau_f)
                tau_f = tau_f + dtau
            enddo
        enddo
    enddo
enddo
end subroutine

subroutine retarded_add(t_lo, t_hi, c_a, c_b, interp)
!add the source function to the images for all cells whose retarded time t - delta_x lies in (t_lo, t_hi]:
!interp = .true.: linear interpolation between s_prev (at time c_a) and s_cur (at time c_b)
!interp = .false.: s_prev
implicit none
real(DP), intent(in) :: t_lo, t_hi, c_a, c_b
logical, intent(in) :: interp
integer :: isel, m, ii, jj, kk, ii_lo, ii_hi
real(DP) :: f(nx)

!$omp parallel do collapse(2) private(ii, jj, kk, ii_lo, ii_hi, f) schedule(dynamic)
do isel=1, n_sel
    do m=1, lc_n
        ii_lo = nx+1; ii_hi = 0
        f = 0
        do ii=1, nx
            if ((delta_x(ii).ge.lc_t(m)-t_hi).and.(delta_x(ii).lt.lc_t(m)-t_lo)) then
                ii_lo = min(ii_lo, ii); ii_hi = max(ii_hi, ii)
                if (interp) f(ii) = min(max((lc_t(m)-delta_x(ii)-c_a)/(c_b-c_a), 0.0_DP), 1.0_DP)
            endif
        enddo
        do kk=1, nz
            do jj=1, ny
                do ii=ii_lo, ii_hi
                    if (w_rt(ii,jj,kk,isel).le.0) cycle
                    img_t(jj,kk,m,isel) = img_t(jj,kk,m,isel) + w_rt(ii,jj,kk,isel) &
                        *((1-f(ii))*s_prev(ii,jj,kk,isel) + f(ii)*s_cur(ii,jj,kk,isel))
                enddo
            enddo
        enddo
    enddo
enddo
!$omp end parallel do
end subroutine

subroutine mrw_init()
!table of the probability to be still inside a sphere of radius R0 for a photon diffusing from its centre:
!P(y) = 2 sum_n (-1)**(n+1) exp(-n**2 pi**2 y), y = D t/R0**2, D = c/(3 rho kappa_ext)
implicit none
integer, parameter :: n_y = 2000
integer :: iy, n
real(DP) :: term
allocate(y_mrw(n_y), p_mrw(n_y))
do iy=1, n_y
    y_mrw(iy) = 1e-4_DP*(5/1e-4_DP)**(real(iy-1,DP)/(n_y-1))
    p_mrw(iy) = 0
    n = 0
    do
        n = n + 1
        term = exp(-n**2*pi**2*y_mrw(iy))
        p_mrw(iy) = p_mrw(iy) + 2*(-1)**(n+1)*term
        if (term.lt.1e-18_DP) exit
    enddo
    p_mrw(iy) = min(max(p_mrw(iy), 0.0_DP), 1.0_DP)
enddo
end subroutine

function mrw_length(chi, r0cm)
!path length [cm] of a diffusing photon until it first reaches the sphere of radius r0cm [cm],
!chi = rho*kappa_ext [cm-1]: draw zeta, solve P(y) = zeta, path = c t = 3 chi R_eff**2 y
implicit none
real(DP), intent(in) :: chi, r0cm
real(DP) :: mrw_length, zeta, y_s, f
integer :: lo, hi, mid
call random_number(zeta)
zeta = max(zeta, tiny(1.0_DP))
if (zeta.le.p_mrw(size(p_mrw))) then
    y_s = -log(zeta/2)/pi**2 !first term of the series, exact enough there
else if (zeta.ge.p_mrw(1)) then
    y_s = y_mrw(1)
else
    lo = 1; hi = size(p_mrw) !p_mrw decreases with y
    do while (hi-lo.gt.1)
        mid = (lo+hi)/2
        if (p_mrw(mid).ge.zeta) then
            lo = mid
        else
            hi = mid
        endif
    enddo
    f = (zeta-p_mrw(lo))/(p_mrw(hi)-p_mrw(lo))
    y_s = y_mrw(lo) + f*(y_mrw(hi)-y_mrw(lo))
endif
!the walk consists of discrete flights: a photon diffusing to radius r0 behaves as if the absorbing sphere
!were larger by the Milne extrapolation length 0.7104 mean free paths (checked with a brute-force walk)
mrw_length = max(3*chi*(r0cm+0.7104_DP/chi)**2*y_s, r0cm) !cannot be shorter than the straight line
end function

function mrw_eff_length(a_abs, l)
!path length weighted with the packet weight, which decays as exp(-a_abs*s) along the walk:
!int_0^l exp(-a_abs s) ds (expected weight after N interactions albedo**N with N ~ Poisson(chi l))
implicit none
real(DP), intent(in) :: a_abs, l
real(DP) :: mrw_eff_length
if (a_abs*l.gt.1e-8_DP) then
    mrw_eff_length = (1-exp(-a_abs*l))/a_abs
else
    mrw_eff_length = l
endif
end function

subroutine random_dir(ux, uy, uz)
!isotropic unit vector
implicit none
real(DP), intent(out) :: ux, uy, uz
real(DP) :: xi1, xi2, phi, cos_theta, sin_theta
call random_number(xi1)
call random_number(xi2)
phi = 2*pi*xi1
cos_theta = 2*xi2 - 1
sin_theta = sqrt(1 - cos_theta**2)
ux = sin_theta*cos(phi)
uy = sin_theta*sin(phi)
uz = cos_theta
end subroutine

subroutine outward_dir(ex, ey, ez, ux, uy, uz)
!direction of a photon crossing a surface with outward normal (ex, ey, ez) in a nearly isotropic
!(diffusive) radiation field: p(mu) ~ mu for 0 < mu < 1, i.e. mu = sqrt(xi)
implicit none
real(DP), intent(in) :: ex, ey, ez
real(DP), intent(out) :: ux, uy, uz
real(DP) :: xi1, xi2, mu, sin_t, phi, ax, ay, az, bx, by, bz, nrm
call random_number(xi1)
call random_number(xi2)
mu = sqrt(xi1)
sin_t = sqrt(1 - mu**2)
phi = 2*pi*xi2
!orthonormal basis (a, b) perpendicular to e
if (abs(ex).lt.0.9) then
    ax = 0; ay = ez; az = -ey     ! e x (1,0,0)
else
    ax = -ez; ay = 0; az = ex     ! e x (0,1,0)
endif
nrm = sqrt(ax**2+ay**2+az**2)
ax = ax/nrm; ay = ay/nrm; az = az/nrm
bx = ey*az - ez*ay; by = ez*ax - ex*az; bz = ex*ay - ey*ax
ux = mu*ex + sin_t*(cos(phi)*ax + sin(phi)*bx)
uy = mu*ey + sin_t*(cos(phi)*ay + sin(phi)*by)
uz = mu*ez + sin_t*(cos(phi)*az + sin(phi)*bz)
end subroutine

elemental function B(nu, T)
implicit none
real(DP), intent(in) :: nu, T ![nu] = Hz, [T] = K, [B] = erg/(cm^2 s Hz sr)
real(DP) :: B, xx
if (T.le.0) then
    B=0
    return
endif
xx = hP*nu/kB/T
if (xx.gt.700) then !avoid overflow of exp, B is negligible there
    B = 0
else
    B = 2*hP*nu**3/clight**2 * 1/(exp(xx)-1)
endif
end function
end program part1_3d
