# Makefile for the 3D Monte Carlo radiative transfer code
#   make          optimised build with OpenMP  -> ./mc3d
#   make debug    build with runtime checks    -> ./mc3d
#   make clean    remove executable and module files

FC     = gfortran
FFLAGS = -O3 -march=native -fopenmp
DFLAGS = -O0 -g -fopenmp -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero

SRC = testmc3d.f90
EXE = mc3d

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(EXE)

debug: $(SRC)
	$(FC) $(DFLAGS) $(SRC) -o $(EXE)

clean:
	rm -f $(EXE) *.mod

.PHONY: debug clean
