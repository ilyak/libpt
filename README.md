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
