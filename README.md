# libpt

This is work-in-progress.

Libpt is a fast hybrid MPI/OpenMP parallel implementation of the CCSD(T) energy
correction for the distributed multiprocessor systems. Due to hybrid design the
code efficiently uses fast memory by sharing data within a single node while
using MPI for inter-node communication. The code exhibits excellent scalability
up to hundreds of cores.

### Code structure

The code was designed to be easily embeddable. Simply include `pt.h` header and
add `pt.c` to the list of project sources. Note that you must compile your
project with MPI and OpenMP enabled.
