ray-sz: ray tracing in Szekeres geometries
=======================================
This is a code to

* Numerically solve the Einstein field equations for the [Szekeres class
of solutions](https://projecteuclid.org/euclid.cmp/1103860587).
This is done by first transforming the equations down to a
small set of elliptic integrals that can be evaluated efficiently as
[Carlson symmetric forms](https://en.wikipedia.org/wiki/Carlson_symmetric_form).
* Ray trace through these spacetimes and return relevant
observational quantities of interest.

The Szekeres solutions represent a general class of exact inhomogeneous
spacetimes. These spacetimes are useful for modelling highly nonlinear
~10 Mpc cosmological structures, such as galaxy clusters or voids, embedded in
a FLRW background.
These solutions encompass the more frequently used spherically symmetric solutions,
inlcuding the LTB family and by extension the FLRW solutions too (though it
is not designed to do so in an efficient way).

A simple Python wrapper is also included. Compile with `make raysz`.


Dependencies
------------
Most of the code is self-contained but we make light use of

* [HEALPix library](http://healpix.sourceforge.net/)

where we need to initalise geodesics on the sky and analyse resulting 
temperature maps. Note in order to install Healpix one will also
need

* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)


