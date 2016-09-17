# gcc with Netlib BLAS on OpenBSD
CC= mpicc
CFLAGS= -std=c99 -Wall -Wextra -g -O3
LDFLAGS= -L/usr/local/lib
LIBS= -lblas -lg2c -lm

# icc with MKL on Linux
#CC= mpicc
#CFLAGS= -std=c99 -Wall -Wextra -g -O3 -mkl=sequential -fopenmp
#LDFLAGS=
#LIBS=

ALL_O= pt.o ptcmd.o strtonum.o
RM= rm -f

all: pt

pt: $(ALL_O)
	$(CC) -o $@ $(CFLAGS) $(ALL_O) $(LDFLAGS) $(LIBS)

check: pt
	@echo rpt01 && ./pt -t tests/rpt01.dat && echo success
	@echo rpt02 && ./pt -t tests/rpt02.dat && echo success
	@echo rpt03 && ./pt -t tests/rpt03.dat && echo success
	@echo rpt04 && ./pt -t tests/rpt04.dat && echo success
	@echo rpt05 && ./pt -t tests/rpt05.dat && echo success
	@echo rpt06 && ./pt -t tests/rpt06.dat && echo success
	@echo rpt07 && ./pt -t tests/rpt07.dat && echo success
	@echo rpt08 && ./pt -t tests/rpt08.dat && echo success
	@echo rpt09 && ./pt -t tests/rpt09.dat && echo success
	@echo upt01 && ./pt -t tests/upt01.dat && echo success
	@echo upt02 && ./pt -t tests/upt02.dat && echo success
	@echo upt03 && ./pt -t tests/upt03.dat && echo success

checkmpi: pt
	@echo rpt01 && mpirun -np 2 ./pt -t tests/rpt01.dat && echo success
	@echo rpt02 && mpirun -np 3 ./pt -t tests/rpt02.dat && echo success
	@echo rpt03 && mpirun -np 4 ./pt -t tests/rpt03.dat && echo success
	@echo rpt04 && mpirun -np 3 ./pt -t tests/rpt04.dat && echo success
	@echo rpt05 && mpirun -np 1 ./pt -t tests/rpt05.dat && echo success
	@echo rpt06 && mpirun -np 2 ./pt -t tests/rpt06.dat && echo success
	@echo rpt07 && mpirun -np 4 ./pt -t tests/rpt07.dat && echo success
	@echo rpt08 && mpirun -np 3 ./pt -t tests/rpt08.dat && echo success
	@echo rpt09 && mpirun -np 2 ./pt -t tests/rpt09.dat && echo success
	@echo upt01 && mpirun -np 3 ./pt -t tests/upt01.dat && echo success
	@echo upt02 && mpirun -np 4 ./pt -t tests/upt02.dat && echo success
	@echo upt03 && mpirun -np 3 ./pt -t tests/upt03.dat && echo success

clean:
	$(RM) $(ALL_O) pt gmon.out *.core *.log

.PHONY: all check clean
