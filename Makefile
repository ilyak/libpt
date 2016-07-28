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
	@./pt -t tests/pt01new.dat && echo success
	@./pt -t tests/pt02new.dat && echo success
	@./pt -t tests/pt03new.dat && echo success
	@./pt -t tests/pt04new.dat && echo success
	@./pt -t tests/pt05new.dat && echo success
	@./pt -t tests/pt06new.dat && echo success
	@./pt -t tests/pt07new.dat && echo success
	@./pt -t tests/pt08new.dat && echo success
	@./pt -t tests/pt09new.dat && echo success

clean:
	$(RM) $(ALL_O) pt gmon.out *.core *.log

.PHONY: all check clean
