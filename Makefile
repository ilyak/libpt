CC= mpicc
CFLAGS= -I../libxutil -Wall -Wextra -g -O3
LDFLAGS= -L../libxutil
#LIBS= -lxutil -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm
LIBS= -lxutil -lblas -lg2c -lm

ALL_O=
RM= rm -f

all: pt

pt: pt.o $(ALL_O)
	$(CC) -o $@ $(CFLAGS) pt.o $(ALL_O) $(LDFLAGS) $(LIBS)

check: pt
	@./pt -t test1.dat && echo success
	@./pt -t test2.dat && echo success
	@./pt -t test3.dat && echo success
	@./pt -t test4.dat && echo success
	@./pt -t test5.dat && echo success

clean:
	$(RM) $(ALL_O) pt pt.o gmon.out *.core *.log

.PHONY: all check clean
