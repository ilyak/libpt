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
	@echo pt01  && ./pt -t tests/pt01.dat && echo success
	@echo pt02  && ./pt -t tests/pt02.dat && echo success
	@echo pt03  && ./pt -t tests/pt03.dat && echo success
	@echo pt04  && ./pt -t tests/pt04.dat && echo success
	@echo pt05  && ./pt -t tests/pt05.dat && echo success
	@echo pt06  && ./pt -t tests/pt06.dat && echo success
	@echo pt07  && ./pt -t tests/pt07.dat && echo success
	@echo pt08  && ./pt -t tests/pt08.dat && echo success
	@echo pt09  && ./pt -t tests/pt09.dat && echo success
	@echo upt01 && ./pt -t tests/upt01.dat && echo success
	@echo upt02 && ./pt -t tests/upt02.dat && echo success
	@echo upt03 && ./pt -t tests/upt03.dat && echo success

checkmpi: pt
	@echo pt01  && mpirun -np 3 ./pt -t tests/pt01.dat && echo success
	@echo pt02  && mpirun -np 3 ./pt -t tests/pt02.dat && echo success
	@echo pt03  && mpirun -np 3 ./pt -t tests/pt03.dat && echo success
	@echo pt04  && mpirun -np 3 ./pt -t tests/pt04.dat && echo success
	@echo pt05  && mpirun -np 3 ./pt -t tests/pt05.dat && echo success
	@echo pt06  && mpirun -np 3 ./pt -t tests/pt06.dat && echo success
	@echo pt07  && mpirun -np 3 ./pt -t tests/pt07.dat && echo success
	@echo pt08  && mpirun -np 3 ./pt -t tests/pt08.dat && echo success
	@echo pt09  && mpirun -np 3 ./pt -t tests/pt09.dat && echo success
	@echo upt01 && mpirun -np 3 ./pt -t tests/upt01.dat && echo success
	@echo upt02 && mpirun -np 3 ./pt -t tests/upt02.dat && echo success
	@echo upt03 && mpirun -np 3 ./pt -t tests/upt03.dat && echo success

clean:
	$(RM) $(ALL_O) pt gmon.out *.core *.log

.PHONY: all check clean
