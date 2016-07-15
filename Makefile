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
	@./pt -t tests/test1.dat && echo success
	@./pt -t tests/test2.dat && echo success
	@./pt -t tests/test3.dat && echo success
	@./pt -t tests/test4.dat && echo success
	@./pt -t tests/test5.dat && echo success
	@./pt -t tests/test6.dat && echo success
	@./pt -t tests/test7.dat && echo success
	@./pt -t tests/test8.dat && echo success
	@./pt -t tests/test9.dat && echo success

clean:
	$(RM) $(ALL_O) pt ptcmd.o gmon.out *.core *.log

.PHONY: all check clean
