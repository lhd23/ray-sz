ray-sz: a general purpose Szekeres code
=======================================
``ray-sz`` is a Fortran code for

* numerically solving the Einstein field equations for the Szekeres family
of solutions by reducing the equations down to a small set of elliptic
integrals that can be evaluated quickly as Carlson R-functions
* raytracing through these spacetimes and return the relevant optical and
observational quantities

Note the Szekeres solutions generalises the LTB family and, by extension,
the FLRW solutions too (though it can recover this special case it is not
designed to do so in an efficient manner).


Dependencies
------------
This code requires the following:

* [HEALPix library](http://healpix.sourceforge.net/)
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)

