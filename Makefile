CC= cc
CFLAGS= -Wall -Wextra -g -fopenmp -Isrc
LDFLAGS= -L/usr/local/lib
LIBS= -lblas -lm

# Intel Compiler with MPI and OpenMP
#CC= mpicc
#CFLAGS= -Wall -Wextra -g -O3 -mkl=sequential -fopenmp -DLIBPT_USE_MPI -Isrc
#LDFLAGS=
#LIBS=

LIBPT_A= src/libpt.a

all: benchmark test

$(LIBPT_A):
	cd src && CC="$(CC)" CFLAGS="$(CFLAGS)" $(MAKE)

benchmark: $(LIBPT_A) benchmark.o
	$(CC) -o $@ $(CFLAGS) benchmark.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

test: $(LIBPT_A) test.o
	$(CC) -o $@ $(CFLAGS) test.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

check: test
	@./test rpt 01 && echo success
	@./test rpt 02 && echo success
	@./test rpt 03 && echo success
	@./test rpt 04 && echo success
	@./test rpt 05 && echo success
	@./test rpt 06 && echo success
	@./test rpt 07 && echo success
	@./test rpt 08 && echo success
	@./test rpt 09 && echo success
	@./test upt 01 && echo success
	@./test upt 02 && echo success
	@./test upt 03 && echo success
	@./test rft 01 && echo success
	@./test rft 02 && echo success
	@./test rft 03 && echo success
	@./test uft 01 && echo success
	@./test uft 02 && echo success
	@./test uft 03 && echo success

checkmpi: test
	@mpirun -np 2 ./test rpt 01 && echo success
	@mpirun -np 3 ./test rpt 02 && echo success
	@mpirun -np 4 ./test rpt 03 && echo success
	@mpirun -np 3 ./test rpt 04 && echo success
	@mpirun -np 1 ./test rpt 05 && echo success
	@mpirun -np 2 ./test rpt 06 && echo success
	@mpirun -np 4 ./test rpt 07 && echo success
	@mpirun -np 3 ./test rpt 08 && echo success
	@mpirun -np 2 ./test rpt 09 && echo success
	@mpirun -np 3 ./test upt 01 && echo success
	@mpirun -np 4 ./test upt 02 && echo success
	@mpirun -np 3 ./test upt 03 && echo success
	@mpirun -np 3 ./test rft 01 && echo success
	@mpirun -np 2 ./test rft 02 && echo success
	@mpirun -np 1 ./test rft 03 && echo success
	@mpirun -np 3 ./test uft 01 && echo success
	@mpirun -np 2 ./test uft 02 && echo success
	@mpirun -np 3 ./test uft 03 && echo success

clean:
	cd src && $(MAKE) clean
	rm -f *.core *.o gmon.out benchmark test

.PHONY: all check checkmpi clean
