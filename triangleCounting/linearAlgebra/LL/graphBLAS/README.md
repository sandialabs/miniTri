# Graph BLAS implementation

This implementation of triangle counting uses the Graph BLAS implementation by Tim Davis.

Currently, this implementation also uses the matrix reading and sorting functionality provided by Kokkos Kernels.
Eventually, we can implement this separately.


## Building Instructions

In order to build this code, you'll need to set the following environment variables
(or change them in the makefile) KOKKOS_PATH and KK_PATH, which should point to the installation of the
Kokkos and Kokkos-Kernels libraries, respectively.  Kokkos and Kokkos-Kernels can be found at
https://github.com/kokkos .

You'll also need to set the Graph BLAS stuff... directions to come.

