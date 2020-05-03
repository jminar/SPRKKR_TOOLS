# GNU Fortran compiler:
#
FC=gfortran
FFLAGS = -cpp
FFLAGS += -fcheck=all -Wall -fbacktrace -fbounds-check -Warray-temporaries -Wconversion-extra

# Intel Fortran compiler:
#
# FG=ifort
# FFLAGS = -fpp
# FFLAGS += -check all -traceback

all: Makefile sum_rules_theory.f
	$(FC) $(FFLAGS) sum_rules_theory.f -o sum_rules_theory.x 

clean:
	rm sum_rules_theory.x *.mod 
