# libpt

Libpt is a fast hybrid MPI/OpenMP parallel implementation of the CCSD(T) energy
correction for the distributed multiprocessor systems. Due to hybrid design the
code efficiently uses fast memory by sharing data within a single node while
using MPI for inter-node communication.

The code exhibits excellent scalability up to hundreds of cores. The benchmark
results are shown in the following table.

Total worker threads | Execution time, min |
---------------------|---------------------|
           8         |         98          |
          32         |         25          |
          64         |         15          |
         128         |          8          |
         256         |          5          |

The benchmark was performed using 1 MPI process per node with 8 threads per
process. The benchmark parameters were `o = 30` and `v = 230`.
