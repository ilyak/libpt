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

all: testrpt testupt testrft testuft benchmarkrpt benchmarkuft

testrpt: pt.o testrpt.o
	$(CC) -o $@ $(CFLAGS) pt.o testrpt.o $(LDFLAGS) $(LIBS)

testupt: pt.o testupt.o
	$(CC) -o $@ $(CFLAGS) pt.o testupt.o $(LDFLAGS) $(LIBS)

testrft: pt.o testrft.o
	$(CC) -o $@ $(CFLAGS) pt.o testrft.o $(LDFLAGS) $(LIBS)

testuft: pt.o testuft.o
	$(CC) -o $@ $(CFLAGS) pt.o testuft.o $(LDFLAGS) $(LIBS)

benchmarkrpt: pt.o benchmarkrpt.o
	$(CC) -o $@ $(CFLAGS) pt.o benchmarkrpt.o $(LDFLAGS) $(LIBS)

benchmarkuft: pt.o benchmarkuft.o
	$(CC) -o $@ $(CFLAGS) pt.o benchmarkuft.o $(LDFLAGS) $(LIBS)

check: testrpt testupt testrft testuft
	@echo rpt01 && ./testrpt tests/rpt01.dat && echo success
	@echo rpt02 && ./testrpt tests/rpt02.dat && echo success
	@echo rpt03 && ./testrpt tests/rpt03.dat && echo success
	@echo rpt04 && ./testrpt tests/rpt04.dat && echo success
	@echo rpt05 && ./testrpt tests/rpt05.dat && echo success
	@echo rpt06 && ./testrpt tests/rpt06.dat && echo success
	@echo rpt07 && ./testrpt tests/rpt07.dat && echo success
	@echo rpt08 && ./testrpt tests/rpt08.dat && echo success
	@echo rpt09 && ./testrpt tests/rpt09.dat && echo success
	@echo upt01 && ./testupt tests/upt01.dat && echo success
	@echo upt02 && ./testupt tests/upt02.dat && echo success
	@echo upt03 && ./testupt tests/upt03.dat && echo success
	@echo rft01 && ./testrft tests/rft01.dat && echo success
	@echo rft02 && ./testrft tests/rft02.dat && echo success
	@echo rft03 && ./testrft tests/rft03.dat && echo success
	@echo uft01 && ./testuft tests/uft01.dat && echo success
	@echo uft02 && ./testuft tests/uft02.dat && echo success
	@echo uft03 && ./testuft tests/uft03.dat && echo success

checkmpi: testrpt testupt testrft testuft
	@echo rpt01 && mpirun -np 2 ./testrpt tests/rpt01.dat && echo success
	@echo rpt02 && mpirun -np 3 ./testrpt tests/rpt02.dat && echo success
	@echo rpt03 && mpirun -np 4 ./testrpt tests/rpt03.dat && echo success
	@echo rpt04 && mpirun -np 3 ./testrpt tests/rpt04.dat && echo success
	@echo rpt05 && mpirun -np 1 ./testrpt tests/rpt05.dat && echo success
	@echo rpt06 && mpirun -np 2 ./testrpt tests/rpt06.dat && echo success
	@echo rpt07 && mpirun -np 4 ./testrpt tests/rpt07.dat && echo success
	@echo rpt08 && mpirun -np 3 ./testrpt tests/rpt08.dat && echo success
	@echo rpt09 && mpirun -np 2 ./testrpt tests/rpt09.dat && echo success
	@echo upt01 && mpirun -np 3 ./testupt tests/upt01.dat && echo success
	@echo upt02 && mpirun -np 4 ./testupt tests/upt02.dat && echo success
	@echo upt03 && mpirun -np 3 ./testupt tests/upt03.dat && echo success
	@echo rft01 && mpirun -np 3 ./testrft tests/rft01.dat && echo success
	@echo rft02 && mpirun -np 2 ./testrft tests/rft02.dat && echo success
	@echo rft03 && mpirun -np 1 ./testrft tests/rft03.dat && echo success
	@echo uft01 && mpirun -np 3 ./testuft tests/uft01.dat && echo success
	@echo uft02 && mpirun -np 2 ./testuft tests/uft02.dat && echo success
	@echo uft03 && mpirun -np 3 ./testuft tests/uft03.dat && echo success

clean:
	rm -f testrpt testupt testrft testuft
	rm -f benchmarkrpt benchmarkuft gmon.out *.core *.o

.PHONY: all check checkmpi clean
