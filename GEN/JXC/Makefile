#
# COMPILER_AND_LIBS := gnu
#
COMPILER_AND_LIBS := intel



ifeq ($(COMPILER_AND_LIBS), intel)
#
# Intel Fortran compiler:
#
  FC=ifort

  FFLAGS = -fpp
  FFLAGS += -O0 -check all -traceback -check all -fpe0 -g -fp-stack-check -ftrapuv -CU
  FFLAGS += -std # Detect integer as logicals etc.;
  FFLAGS += -warn unused # -diag-error 7712 # Really force error on unused variables. The number to trigger this comes from the warning you're trying to elevate to an error.
  FFLAGS += -diag-error 6843 # Really force error on unassigned 'intent(OUT)' variables. The number to trigger this comes from the warning you're trying to elevate to an error.
  # FFLAGS += -diag-once-6916 # Display only once:  warning #6916: Fortran 2008 does not allow this length specification.
else
#
# GNU Fortran compiler:
#
  FC=gfortran
  FFLAGS = -cpp
  FFLAGS += -fcheck=all -Wall -fbacktrace -fbounds-check -Warray-temporaries -Wconversion-extra
# FFLAGS += -pedantic-errors
endif



all: recreate_mag_cluster.F90 Makefile
	$(FC) $(FFLAGS) recreate_mag_cluster.F90 -o recreate_mag_cluster.x

clean:
	rm recreate_mag_cluster.x
