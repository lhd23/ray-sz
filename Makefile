F90C = gfortran
F90FLAGS = -fopenmp -O3 -fimplicit-none -fbacktrace #-ffpe-trap=zero,overflow,underflow #-Wall -g #-no-ipo -ffixed-line-length-none

HEALPIX = /usr/local/Healpix_3.30
CFITSIO = /usr/local/lib

LIBS = -L$(HEALPIX)/lib-gfortran  -L$(CFITSIO) -lhealpix -lcfitsio
INCLUDES = -I$(HEALPIX)/include-gfortran

OBJ = precision1.o constants.o abscissae.o utils.o cosmo_params.o \
		szclass.o szfuncs.o szelliptic.o szlocal.o szgeom.o \
		sznullgeo.o test_module.o


main: $(OBJ) raytrace.o
	$(F90C) $^ $(LIBS) $(INCLUDES) $(F90FLAGS) -o $@

test: $(OBJ) test.o
	$(F90C) $^ $(LIBS) $(INCLUDES) $(F90FLAGS) -o $@

%.o: %.f90
	$(F90C) $(LIBS) $(INCLUDES) $(F90FLAGS) -c $*.f90

test_carlson: $(OBJ) test_carlson.o
	$(F90C) $^ $(F90FLAGS) -o $@


clean: 
	rm -f *.o *.a *.mod test main
