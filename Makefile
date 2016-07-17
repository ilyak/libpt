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
	@./pt -t tests/test01.dat && echo success
	@./pt -t tests/test02.dat && echo success
	@./pt -t tests/test03.dat && echo success
	@./pt -t tests/test04.dat && echo success
	@./pt -t tests/test05.dat && echo success
	@./pt -t tests/test06.dat && echo success
	@./pt -t tests/test07.dat && echo success
	@./pt -t tests/test08.dat && echo success
	@./pt -t tests/test09.dat && echo success

clean:
	$(RM) $(ALL_O) pt ptcmd.o gmon.out *.core *.log

.PHONY: all check clean
