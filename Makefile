F90C = gfortran
F90FLAGS = -fopenmp -O3

# Warning flags
#F90FLAGS += -Wall -Wextra -pedantic -fimplicit-none

# Debug flags
F90FLAGS += -g -fbacktrace -ffpe-trap=zero,overflow,underflow

# Other flags
#F90FLAGS += -ffixed-line-length-none

HEALPIX = /usr/local/Healpix_3.30
CFITSIO = /usr/local/lib

LIBS = -L$(HEALPIX)/lib-gfortran  -L$(CFITSIO) -lhealpix -lcfitsio
INCLUDES = -I$(HEALPIX)/include-gfortran
OMPLIB = -lgomp

LIBDIR = lib
INCDIR = include

SRC = precision1.f90 constants.f90 abscissae.f90 utils.f90 cosmo_params.f90 \
		szclass.f90 szfuncs.f90 elliptic_cache.f90 szelliptic.f90 \
		szlocal.f90 szgeom.f90 sznullgeo.f90 test_module.f90 szcmb.f90
OBJ = $(SRC:.f90=.o)
MOD = $(SRC:.f90=.mod)

all: $(LIBDIR)/libmylib.a main

$(LIBDIR)/libmylib.a: $(SRC)
	mkdir -p $(INCDIR) $(LIBDIR)
	$(F90C) $(LIBS) $(INCLUDES) $(F90FLAGS) -c $(SRC)
	ar -rcs $@ $(OBJ)
	mv $(MOD) $(INCDIR)/
	rm -f $(OBJ)

INCLUDES += -I$(INCDIR)
LIBS += -L$(LIBDIR) -lmylib

main: driver.o
	$(F90C) $^ $(LIBS) $(INCLUDES) $(F90FLAGS) -o $@

raysz:
	f2py -c --fcompiler=gnu95 --f90flags="-fopenmp -O3" $(LIBS) $(OMPLIB) $(INCLUDES) -m $@ raysz.f90

test: $(OBJ) test.o
	$(F90C) $^ $(LIBS) $(INCLUDES) $(F90FLAGS) -o $@

%.o: %.f90
	$(F90C) $(LIBS) $(INCLUDES) $(F90FLAGS) -c $*.f90

clean:
	rm -f *.o *.a *.mod main test $(LIBDIR)/libmylib.a
