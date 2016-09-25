# libpt

Libpt is a fast hybrid MPI/OpenMP parallel implementation of the
coupled-cluster (T) energy correction for clusters of computers and SMP
systems. Due to hybrid design the code efficiently uses RAM by sharing data
within a single node while using MPI only for inter-node communication.

### Code structure

The code was designed to be easily embeddable. Simply include `pt.h` header and
add `pt.c` to the list of project sources. OpenMP is enabled automatically when
the appropriate OpenMP compiler flags are used. Defining `WITH_MPI` during
compilation enables MPI support. See `pt.h` file for detailed API
documentation.
