# Kokkos-Kernels implementation

This implementation of triangle counting uses Kokkos and Kokkos-Kernels to count the
number of triangles in parallel.  Currently this works with OpenMP but could be extended 
to work on GPUs as well.

This implementation was similar to implementation that was awarded champion status in
the 2017 IEEE HPEC/DARPA/Amazon Graph Challenge:

> Wolf, Deveci, Berry, Hammond, Rajamanickam. “Fast Linear Algebra-Based Triangle Counting with KokkosKernels,” Proc of IEEE HPEC, 2017.

## Building Instructions

In order to build this code, you'll need to set the following environment variables
(or change them in the makefile) KOKKOS_PATH and KK_PATH, which should point to the installation of the
Kokkos and Kokkos-Kernels libraries, respectively.  Kokkos and Kokkos-Kernels can be found at
https://github.com/kokkos .


