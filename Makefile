# gcc with Netlib BLAS on OpenBSD
CC= cc
CFLAGS= -Wall -Wextra -g -O3
LDFLAGS= -L/usr/local/lib
LIBS= -lblas -lg2c -lm

# gcc with Netlib BLAS and MPI on OpenBSD
#CC= mpicc
#CFLAGS= -Wall -Wextra -g -O3 -DWITH_MPI
#LDFLAGS= -L/usr/local/lib
#LIBS= -lblas -lg2c -lm

# icc with MKL and MPI on Linux
#CC= mpicc
#CFLAGS= -Wall -Wextra -g -O3 -mkl=sequential -fopenmp -DWITH_MPI
#LDFLAGS=
#LIBS=

all: testpt testft

testpt: pt.o testpt.o
	$(CC) -o $@ $(CFLAGS) pt.o testpt.o $(LDFLAGS) $(LIBS)

testft: pt.o testft.o
	$(CC) -o $@ $(CFLAGS) pt.o testft.o $(LDFLAGS) $(LIBS)

check: testpt testft
	@echo rpt01 && ./testpt -t tests/rpt01.dat && echo success
	@echo rpt02 && ./testpt -t tests/rpt02.dat && echo success
	@echo rpt03 && ./testpt -t tests/rpt03.dat && echo success
	@echo rpt04 && ./testpt -t tests/rpt04.dat && echo success
	@echo rpt05 && ./testpt -t tests/rpt05.dat && echo success
	@echo rpt06 && ./testpt -t tests/rpt06.dat && echo success
	@echo rpt07 && ./testpt -t tests/rpt07.dat && echo success
	@echo rpt08 && ./testpt -t tests/rpt08.dat && echo success
	@echo rpt09 && ./testpt -t tests/rpt09.dat && echo success
	@echo upt01 && ./testpt -t tests/upt01.dat && echo success
	@echo upt02 && ./testpt -t tests/upt02.dat && echo success
	@echo upt03 && ./testpt -t tests/upt03.dat && echo success

checkmpi: testpt testft
	@echo rpt01 && mpirun -np 2 ./testpt -t tests/rpt01.dat && echo success
	@echo rpt02 && mpirun -np 3 ./testpt -t tests/rpt02.dat && echo success
	@echo rpt03 && mpirun -np 4 ./testpt -t tests/rpt03.dat && echo success
	@echo rpt04 && mpirun -np 3 ./testpt -t tests/rpt04.dat && echo success
	@echo rpt05 && mpirun -np 1 ./testpt -t tests/rpt05.dat && echo success
	@echo rpt06 && mpirun -np 2 ./testpt -t tests/rpt06.dat && echo success
	@echo rpt07 && mpirun -np 4 ./testpt -t tests/rpt07.dat && echo success
	@echo rpt08 && mpirun -np 3 ./testpt -t tests/rpt08.dat && echo success
	@echo rpt09 && mpirun -np 2 ./testpt -t tests/rpt09.dat && echo success
	@echo upt01 && mpirun -np 3 ./testpt -t tests/upt01.dat && echo success
	@echo upt02 && mpirun -np 4 ./testpt -t tests/upt02.dat && echo success
	@echo upt03 && mpirun -np 3 ./testpt -t tests/upt03.dat && echo success

clean:
	rm -f testpt testft gmon.out *.core *.log *.o

.PHONY: all check checkmpi clean
