# libpt

Libpt is a fast MPI/OpenMP parallel implementation of the coupled-cluster (T)
and (fT) energy corrections for clusters of computers and SMP systems. Due to
hybrid parallel design the code efficiently uses RAM by sharing data within a
single node while using message-passing only for inter-node communication.

The API defines several functions:

- `libpt_rpt` - computes (T) energy correction for the restricted case
- `libpt_upt` - computes (T) energy correction for the unrestricted case
- `libpt_rft` - computes (fT) energy correction for the restricted case
- `libpt_uft` - computes (fT) energy correction for the unrestricted case

All functions are MPI/OpenMP parallel. The detailed documentation can be found
in the `pt.h` file. To use libpt in your project, include `pt.h` header and add
`pt.c` to the list of project source files. OpenMP is enabled automatically
when appropriate OpenMP compiler flags are used. Defining `WITH_MPI` during
compilation enables MPI support.

### Benchmarks

Libpt shows excellent performance both on single workstations as well as on
distributed clusters of computers. Parallel scalability benchmarks on the NERSC
Cori Cray XC40 supercomputer (32 CPU cores per node) are shown below. A cluster
of 9 water molecules with the cc-pVTZ basis set (522 basis functions in total)
was used as a benchmark system.

|      Nodes       | (T) calculation time | Speedup (max) |
|:----------------:|:--------------------:|:-------------:|
|   4  (128 cores) |            11064 sec |      1x (1x)  |
|   8  (256 cores) |             5599 sec |      2x (2x)  |
|  16  (512 cores) |             2925 sec |    3.8x (4x)  |
|  32 (1024 cores) |             1487 sec |    7.4x (8x)  |
|  64 (2048 cores) |              819 sec |   13.5x (16x) |
| 128 (4096 cores) |              440 sec |   25.1x (32x) |
