CC= mpicc
CFLAGS= -Wall -Wextra -g -O3
LDFLAGS=
LIBS= -lblas -lg2c -lm

#CC= mpicc
#CFLAGS= -Wall -Wextra -g -O3 -fopenmp
#LDFLAGS=
#LIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

ALL_O= pt.o ptcmd.o strtonum.o reallocarray.o
RM= rm -f

all: pt

pt: $(ALL_O)
	$(CC) -o $@ $(CFLAGS) $(ALL_O) $(LDFLAGS) $(LIBS)

check: pt
	@./pt -t tests/pt01.dat && echo success
	@./pt -t tests/pt02.dat && echo success
	@./pt -t tests/pt03.dat && echo success
	@./pt -t tests/pt04.dat && echo success
	@./pt -t tests/pt05.dat && echo success
	@./pt -t tests/pt06.dat && echo success
	@./pt -t tests/pt07.dat && echo success
	@./pt -t tests/pt08.dat && echo success
	@./pt -t tests/pt09.dat && echo success

clean:
	$(RM) $(ALL_O) pt gmon.out *.core *.log

.PHONY: all check clean
