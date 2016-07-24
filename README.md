# libpt

This is work-in-progress.

Libpt is a fast hybrid MPI/OpenMP parallel implementation of the CCSD(T) energy
correction for the distributed multiprocessor systems. Due to hybrid design the
code efficiently uses fast memory by sharing data within a single node while
using MPI for inter-node communication. The code exhibits excellent scalability
up to hundreds of cores.

Total worker threads | Execution time, min |
---------------------|---------------------|
           8         |         98          |
          32         |         25          |
          64         |         15          |
         128         |          8          |
         256         |          5          |

The benchmark shown above was performed on a computer cluster with Ethernet
interconnect. Each node had two Quadcore AMD Opteron 2.3 GHz processors with
16 Gb of main memory. One MPI process per node with 8 threads per process were
used. Benchmark size parameters were `o = 30` and `v = 230`.

### Code structure

The code was designed to be easily embeddable. Simply include `pt.h` header and
add `pt.c` to the list of project sources. Note that you must compile your
project with MPI and OpenMP enabled.
