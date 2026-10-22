# testmc

`testmc` started as a simple study project for my theoretical astrophysics course (in 2022). It is a gross simplification of a [RADMC-3D](https://github.com/dullemond/radmc3d-2.0) functionality written following [Dullemond's lectures](https://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2012/) for method demonstration purposes and now it is repurposed for testing Monte Carlo RT with various addons (like time-dependent light propagation). The original version is `testmc2d.f90` and `testmc3d.f90` is a vibecoded expansion.

## testmc2d.f90

A ready to go 2D Monte Carlo RT for making monochromatic isotropically scattered light images. It is hardcoded with a spherical cloud in Cartesian coordinates, a star in the centre and three wavelengths and Qsca that I got I don't remember where. There are two options to compute the image: either first the source function is estimated and then it is ray-traced, or you can make an MC image directly. You simply compile it with your favourite Fortran compiler and enjoy your 1D intensity images and 2D source function at all 3 wavelengths. All changes must be made inside .f90.

## testmc3d.f90

Computah, make my code 3D. And it did. Additionally, it reads an opacity file of the same structure like RADMC-3D (only one dust species, scattering only isotropic and always named `dustkappa.inp` in the folder where you run the code) and is parallelised with OpenMP (takes all threads it can get but you can limit it through environment variables). Absorption is now included so energy is removed from packets and they are handled using survival biasing and Russian roulette methods. The images is still just scattered light, however. The output of source function was replaced by temperature in a dedicated mode. The dust is, however, not self-heated so the temperature will be underestimated in dense regions. The code can be compiled by running `make`. `make debug` will compile with runtime checks and `make clean` will remove the executable. Once compiled it can be called with:

```bash
./mc3d imagert inu 3 2 5 nphot 1e4
```

Possible modes are: 

`imagert` --- MC source function + ray-tracing (default mode), output: image_raytrace_inu.out, where inu are wavelength indices of the grid in dustkappa.inp

Inside is a 2D array of flux/pixel in erg s-1 Hz-1 cm-2 with observer at -150 pc on x-axis, rows are along z-coordinate and columns are along y-coordinate. 

`imagemc` --- image directly from MC, output is the same

`mctherm` --- MC source function but output is dust_temperature3d.out with 3D array of temperature in Kelvin and is written in XY-blocks, rows are along y-coordinate, columns are along x-coordinate and the matrix repeats nz times.

There are two additional keywords to be used in the terminal: `inu` restricts which wavelengths to output (by default it outputs all) and `nphot` sets the amount of photon packages per wavelength(!).

The hardcoded grid (Cartesian 101x101x101 with 1 au step) and density (spherically symmetric cloud with rho=rho0(1au/r)² from 1 to 40 au) remain and if you wish a more refined grid or different density, you have to change .f90 directly and recompile.
