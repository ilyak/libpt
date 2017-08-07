CC= cc
CFLAGS= -Wall -Wextra -g -fopenmp -Isrc
LDFLAGS= -L/usr/local/lib
LIBS= -lblas -lm

# Intel Compiler with MPI and OpenMP
#CC= mpicc
#CFLAGS= -Wall -Wextra -g -O3 -mkl=sequential -fopenmp -DLIBPT_USE_MPI -Isrc
#LDFLAGS=
#LIBS=

ALL_TESTS= testrpt testupt testrft testuft
ALL_BENCHMARKS= benchmarkrpt benchmarkupt benchmarkrft benchmarkuft

LIBPT_A= src/libpt.a

all: $(ALL_TESTS) $(ALL_BENCHMARKS)

$(LIBPT_A):
	cd src && CC="$(CC)" CFLAGS="$(CFLAGS)" $(MAKE)

testrpt: $(LIBPT_A) testpt.o
	$(CC) -o $@ $(CFLAGS) testpt.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

testupt: $(LIBPT_A) testpt.o
	$(CC) -o $@ $(CFLAGS) testpt.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

testrft: $(LIBPT_A) testft.o
	$(CC) -o $@ $(CFLAGS) testft.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

testuft: $(LIBPT_A) testft.o
	$(CC) -o $@ $(CFLAGS) testft.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

benchmarkrpt: $(LIBPT_A) benchmarkpt.o
	$(CC) -o $@ $(CFLAGS) benchmarkpt.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

benchmarkupt: $(LIBPT_A) benchmarkpt.o
	$(CC) -o $@ $(CFLAGS) benchmarkpt.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

benchmarkrft: $(LIBPT_A) benchmarkft.o
	$(CC) -o $@ $(CFLAGS) benchmarkft.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

benchmarkuft: $(LIBPT_A) benchmarkft.o
	$(CC) -o $@ $(CFLAGS) benchmarkft.o $(LDFLAGS) $(LIBPT_A) $(LIBS)

check: $(ALL_TESTS)
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

checkmpi: $(ALL_TESTS)
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
	cd src && $(MAKE) clean
	rm -f *.core *.o gmon.out
	rm -f $(ALL_TESTS) $(ALL_BENCHMARKS)

.PHONY: all check checkmpi clean
