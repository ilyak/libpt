CC= mpicc
CFLAGS= -I../libxutil -Wall -Wextra -g -O3
LDFLAGS= -L../libxutil
LIBS= -lxutil -lblas -lg2c -lm

#CC= mpicc
#CFLAGS= -I../libxutil -Wall -Wextra -g -O3 -fopenmp
#LDFLAGS= -L../libxutil
#LIBS= -lxutil -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

ALL_O= pt.o
RM= rm -f

all: pt

pt: ptcmd.o $(ALL_O)
	$(CC) -o $@ $(CFLAGS) ptcmd.o $(ALL_O) $(LDFLAGS) $(LIBS)

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
	$(RM) $(ALL_O) pt ptcmd.o gmon.out *.core *.log

.PHONY: all check clean
