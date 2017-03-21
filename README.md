# libpt

Libpt is a fast MPI/OpenMP parallel implementation of the coupled-cluster (T)
and (fT) energy corrections for clusters of computers and SMP systems. Due to
hybrid parallel design the code efficiently uses RAM by sharing data within a
single node while using message-passing only for inter-node communication.

To use libpt in your project, include `pt.h` header and add `pt.c` to the list
of project source files. OpenMP is enabled automatically when appropriate
OpenMP compiler flags are used. Defining `WITH_MPI` during compilation enables
MPI support. Detailed API documentation can be found in `pt.h` file.
