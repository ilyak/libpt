# libpt

Libpt is a fast MPI/OpenMP parallel implementation of the coupled-cluster
CCSD(T) and EOM-CCSD(fT) energy corrections for clusters of computers and SMP
systems. Libpt is the default implementation of perturbative triples
corrections in the [Q-Chem](http://www.q-chem.com) package.

### Usage

The API defines several functions:

- `libpt_rpt` - compute CCSD(T) energy correction for the restricted case
- `libpt_upt` - compute CCSD(T) energy correction for the unrestricted case
- `libpt_rft` - compute EOM-CCSD(fT) energy correction for the restricted case
- `libpt_uft` - compute EOM-CCSD(fT) energy correction for the unrestricted case
- `libpt_rpt_mp` - mixed-precision variant of the `libpt_rpt` function
- `libpt_upt_mp` - mixed-precision variant of the `libpt_upt` function
- `libpt_rft_mp` - mixed-precision variant of the `libpt_rft` function
- `libpt_uft_mp` - mixed-precision variant of the `libpt_uft` function

### Compilation

You will need a POSIX-compatible `make` utility, an efficient BLAS library (for
linking), and an ANSI C complaint compiler. To compile libpt, issue:

    cd src && make

To change the default compiler and enable OpenMP and MPI, the following command
can be used:

    cd src && CC=mpicc CFLAGS="-O3 -fopenmp -DLIBPT_USE_MPI" make

To use libpt in your project, include the `pt.h` file and link with the
compiled static `libpt.a` library. The detailed API documentation can be found
in the `pt.h` file.

### Benchmarks

Libpt shows excellent performance both on single workstations as well as on
distributed clusters of computers. Parallel scalability benchmarks on the NERSC
Cori Cray XC40 supercomputer (32 CPU cores per node) are shown below. A cluster
of 9 water molecules with the cc-pVTZ basis set (522 basis functions in total)
was used as a benchmark system.

|      Nodes       | (T) calculation time | Speedup (max) |
|:----------------:|:--------------------:|:-------------:|
|   4  (128 cores) |            11064 sec |      1x  (1x) |
|   8  (256 cores) |             5599 sec |      2x  (2x) |
|  16  (512 cores) |             2925 sec |    3.8x  (4x) |
|  32 (1024 cores) |             1487 sec |    7.4x  (8x) |
|  64 (2048 cores) |              819 sec |   13.5x (16x) |
| 128 (4096 cores) |              440 sec |   25.1x (32x) |
