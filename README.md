# testmc

`testmc` started as a simple study project for my theoretical astrophysics course (in 2022). It is a gross simplification of a [RADMC-3D](https://github.com/dullemond/radmc3d-2.0) functionality written following [Dullemond's lectures](https://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2012/) for method demonstration purposes and now it is repurposed for testing Monte Carlo RT with various addons (like time-dependent light propagation). The original version is `testmc2d.f90` and `testmc3d.f90` is a vibecoded expansion.

## testmc2d.f90

A ready to go 2D Monte Carlo RT for making monochromatic isotropically scattered light images. It is hardcoded with a spherical cloud in Cartesian coordinates, a star in the centre and three wavelengths and Qsca that I got I don't remember where. There are two options to compute the image: either first the source function is estimated and then it is ray-traced, or you can make an MC image directly. You simply compile it with your favourite Fortran compiler and enjoy your 1D intensity images and 2D source function at all 3 wavelengths. All changes must be made inside .f90.

## testmc3d.f90

Computah, make my code 3D. And it did. Additionally, it reads an opacity file of the same structure like RADMC-3D (only one dust species, scattering only isotropic and always named `dustkappa.inp` in the folder where you run the code) and is parallelised with OpenMP (takes all threads it can get but you can limit it through environment variables, e.g. `OMP_NUM_THREADS=8`). Absorption is included: the weight of a packet is reduced by the albedo at every interaction (survival biasing) and packets with very low weight are killed or boosted by Russian roulette. The images can be just scattered light or scattered light+dust thermal emission depending on the mode. The dust is heated by the star (direct and scattered light) and by the CMB (2.73 K), but not by its own emission, so the temperature will be underestimated in regions that are optically thick in the infrared. In optically thick cells the scatterings are replaced by a modified random walk (Min et al. 2009), see below.

The code can be compiled by running `make`. `make debug` will compile with runtime checks and `make clean` will remove the executable. Once compiled it can be called with:

```bash
./mc3d imagert inu 3 2 5 nphot 1e4
```

`./mc3d -h` prints the usage.

### Modes

`imagert` --- MC scattering source function + ray-tracing (default mode), output: `image_raytrace_*.out`, where * are wavelength indices of the grid in `dustkappa.inp`.

`imagemc` --- image directly from MC (photon packets escaping towards the observer within a narrow cone), output is called `image_montecarlo_*.out`. It is noisier than `imagert`.

`mctherm` --- dust temperature, output is `dust_temperature3d.out` in K.

`imagetherm` --- dust temperature + ray-tracing of scattering and dust thermal emission. Output is both `dust_temperature3d.out` and `image_raytrace_therm_*.out`. The image wavelengths have to be given with `inu`.

There are two additional keywords to be used in the terminal: `inu` selects the wavelengths that are computed and written out, and `nphot` sets the number of photon packets. Keywords can come in any order.

* In `imagert` and `imagemc`, every selected wavelength gets `nphot` packets (default: all wavelengths).
* In `mctherm` and `imagetherm`, the dust is heated by `nphot` heating packets in total, whose wavelengths are drawn from the stellar spectrum so that they all carry the same energy. They add up the absorbed power per cell directly, so the cost and memory do not grow with the number of wavelengths in `dustkappa.inp`. Every selected wavelength (`inu`) gets another `nphot` packets for its J and images. Without `inu`, `mctherm` only computes temperatures.
* In time-dependent runs these numbers are per time step.

At the start the code prints the memory needed by the large arrays (and on Linux compares it with the available memory). In time-dependent runs the memory of the packets in flight is printed at every step, and the run stops with an explanation if they would not fit into the available memory. Long loops print their progress.

### Input files

All input files are read from the folder where the code is run. Lines starting with `#` or `!` are comments.

`dustkappa.inp` (required) --- dust opacities in the RADMC-3D format: `iformat` (1: lambda, kappa_abs; 2: + kappa_sca; 3: + g, which is read but ignored), number of wavelengths, then one line per wavelength with lambda in micron and kappa in cm2 per gram of dust. Wavelength indices (`inu`) refer to the lines of this file. For `mctherm`/`imagetherm` the wavelengths must be monotonic and cover both the stellar spectrum and the dust emission (e.g. 0.1-1000 micron).

`grid.inp` and `dust_density.inp` (optional, always together) --- the model. If neither is present, the default model is used: a Cartesian grid of 101x101x101 cells of 1 au and a spherically symmetric shell with rho = rho0 (1 au/r)² from 1 to 40 au. If only one of them is present the code stops and explains both formats.

* `grid.inp`: three rows with the cell edges in au along x, y and z. The edges must increase strictly; the cells do not have to be equal and the three axes can have different numbers of cells (nx = number of x edges - 1, etc.). The star (at x_star, y_star, z_star = 0 in the code) must lie inside the grid.
* `dust_density.inp`: dust mass density in g cm-3 in the same layout as `dust_temperature3d.out` (see below), i.e. ny*nz lines with nx values each. From a numpy array `rho[x, y, z]` it is written with `np.savetxt("dust_density.inp", rho[:, ::-1, :].transpose(2, 1, 0).reshape(-1, rho.shape[0]))`.

`make_disk.py` is an example that writes both files for a disk with its midplane at x = 0 (seen face-on), with fine vertical cells to resolve the scale height.

`lightcurve.inp` (optional) --- if present, the run is time-dependent, see below.

### Output files

Images (`image_*.out`) --- 2D arrays of flux per pixel in erg s-1 cm-2 Hz-1 for an observer at 150 pc on the -x axis. A pixel is one (y, z) cell of the grid. Each line has constant z (the first line is the highest z), the values along a line go from low to high y.

3D arrays (`dust_temperature3d.out`, the J dumps of time-dependent runs, and `dust_density.inp`) --- nz blocks, one per z cell from low to high z. Each block has ny lines, from the HIGHEST y down to the lowest, and every line holds nx values with constant y and z, from low to high x. In numpy:

```python
T = np.loadtxt('dust_temperature3d.out').reshape(nz, ny, nx)[:, ::-1, :].transpose(2, 1, 0)   # T[x, y, z]
img = np.loadtxt('image_raytrace_0003.out')[::-1, :].T                                         # img[y, z]
```

### Time-dependent runs

If `lightcurve.inp` is present, every mode runs time-dependent. Its structure:

```
iformat        1: L_star in Lsun (spectral shape of T_star from the code), 2: R_star in Rsun and T_star in K
nt             number of rows below
itime  time[s]  L_star            (iformat = 1)
itime  time[s]  R_star  T_star    (iformat = 2)
```

The stellar luminosity is interpolated linearly between the rows. If `itime` starts with 0, the run starts from the equilibrium radiation field for the first row. If it starts with 1, the run starts at time 0 from a dark cloud with the dust at 2.73 K, and the first time has to be > 0. Photon packets carry a time stamp. Between consecutive rows they are propagated with a census of the packets in flight, and new packets are emitted during every step. J averaged over each step is dumped to `J_timeseries/J_itime.out`. Each file holds the selected wavelengths (`inu`) one after another, each as a full 3D array in the layout above. In `mctherm` and `imagetherm` the temperature of every step goes to `T_timeseries/T_itime.out`.

Images are made for every row: `image_raytrace_inu_itime.out`, `image_montecarlo_inu_itime.out` and `image_raytrace_therm_inu_itime.out`. Observer time is defined so that the direct starlight has zero delay. `imagert` and `imagetherm` give the image at that instant by ray tracing in retarded time (each cell is seen at t - x/c, with the step averages interpolated linearly between the middles of the steps). `imagemc` gives the flux averaged over the interval since the previous row. After the last row the run continues, with the light curve kept flat, for one light crossing time of the grid, so that the last images are complete.

The time resolution is set by the steps. They should be short compared with the light travel time across the structures you care about (1 au/c is about 500 s). Only the selected wavelengths get packets in `imagert`/`imagemc`; the thermal modes need all wavelengths. The J dumps are large (about 12 MB per wavelength and step for 101³ cells), so use `inu` for long runs. With an equilibrium start, the equilibrium run uses up to 20 times the packets of a time step when the first time step is short compared with the light crossing time of the grid.

### Optically thick cells

Photons in optically thick cells would scatter for a very long time. Where the sphere around a packet that reaches the nearest cell wall has an optical depth above `mrw_gamma` (10 in the code), all scatterings inside this sphere are replaced by a single diffusion step (modified random walk, Min et al. 2009). The path length is drawn from the diffusion solution, with the sphere enlarged by the Milne extrapolation length of 0.71 mean free paths, which was needed to match a brute-force random walk. The packet is then placed on the sphere and leaves it outwards. With absorption, the weight drops by exp(-rho kappa_abs l) over the walk. In time-dependent runs a walk can last several time steps. The walk makes optically thick models feasible, but it does not resolve the structure inside a cell: use a fine enough grid where it matters (e.g. the vertical structure of a disk). In time-dependent runs with an equilibrium start, radiation that stays trapped in very thick, (almost) non-absorbing cells for longer than the age range of the equilibrium packets (1.5 light crossing times of the grid) is part of J at itime 0 but not of the packets in flight. J inside such cells then drops after itime 0 and builds up again on the diffusion time. The code prints a note when this happens. Images and optically thin regions are not affected.
