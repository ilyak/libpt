CFLAGS= -Wall -Wextra -g -fopenmp -Isrc
CXXFLAGS= $(CFLAGS)
LDFLAGS= -L/usr/local/lib
LIBS= -lblas -lm

# Intel Compiler with MPI and OpenMP
#CC= mpicc
#CFLAGS= -Wall -Wextra -g -O3 -mkl=sequential -fopenmp -DLIBPT_USE_MPI -Isrc
#CXXFLAGS= $(CFLAGS)
#LDFLAGS=
#LIBS=

LIBPT_A= src/libpt.a

all: benchmark test

$(LIBPT_A):
	cd src && CC="$(CC)" CFLAGS="$(CFLAGS)" $(MAKE)

benchmark: $(LIBPT_A) benchmark.o
	$(CC) -o $@ $(CFLAGS) benchmark.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

test: $(LIBPT_A) test.o
	$(CXX) -o $@ $(CFLAGS) test.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

check: test
	@./test

checkmpi: test
	@mpirun -np 3 ./test
	@mpirun -np 2 ./test
	@mpirun -np 1 ./test

clean:
	cd src && $(MAKE) clean
	rm -f *.core *.o gmon.out benchmark test

.PHONY: all check checkmpi clean
