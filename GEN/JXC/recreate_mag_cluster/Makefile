# GNU Fortran compiler:
#
# FC=gfortran
# FFLAGS = -cpp
# FFLAGS += -fcheck=all -Wall -fbacktrace -fbounds-check -Warray-temporaries -Wconversion-extra

# Intel Fortran compiler:
#
FC=ifort
FFLAGS = -fpp
FFLAGS += -check all -traceback

all: recreate_mag_cluster.F90 Makefile
	$(FC) $(FFLAGS) recreate_mag_cluster.F90 -o recreate_mag_cluster.x

clean:
	rm recreate_mag_cluster.x
