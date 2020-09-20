ray-sz: ray tracing through exact inhomogeneous spacetimes
=======================================
``ray-sz`` can be used to

* Numerically solve the Einstein field equations for the Szekeres class
of solutions. This is done by first transforming the equations down to a
small set of elliptic integrals that can be evaluated efficiently as
[Carlson symmetric forms](https://en.wikipedia.org/wiki/Carlson_symmetric_form)
* Ray trace through these spacetimes and return the relevant optical and
observational quantities

The Szekeres class of solutions encompasses spherically symmetric solutions,
inlcuding the LTB family and, by extension, the FLRW solutions too (though it
can recover this special case it is not designed to do so in an efficient way).


Dependencies
------------
Most of the code is self-contained but we make light use of

* [HEALPix library](http://healpix.sourceforge.net/)

where we need to initalise geodesics on the sky and compute 
harmonic statistics. Note in order to install Healpix one will also
need

* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)

