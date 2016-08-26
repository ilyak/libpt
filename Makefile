CC= cc
CFLAGS= -std=c99 -Wall -Wextra -g -O3
LDFLAGS= -L/usr/local/lib
LIBS= -lblas -lg2c -lm

#CC= icc
#CFLAGS= -std=c99 -Wall -Wextra -g -O3 -mkl=sequential -fopenmp
#LDFLAGS=
#LIBS=

ALL_O= pt.o ptcmd.o strtonum.o
RM= rm -f

all: pt

pt: $(ALL_O)
	$(CC) -o $@ $(CFLAGS) $(ALL_O) $(LDFLAGS) $(LIBS)

check: pt
	@echo pt01 && ./pt -t tests/pt01.dat && echo success
	@echo pt02 && ./pt -t tests/pt02.dat && echo success
	@echo pt03 && ./pt -t tests/pt03.dat && echo success
	@echo pt04 && ./pt -t tests/pt04.dat && echo success
	@echo pt05 && ./pt -t tests/pt05.dat && echo success
	@echo pt06 && ./pt -t tests/pt06.dat && echo success
	@echo pt07 && ./pt -t tests/pt07.dat && echo success
	@echo pt08 && ./pt -t tests/pt08.dat && echo success
	@echo pt09 && ./pt -t tests/pt09.dat && echo success

clean:
	$(RM) $(ALL_O) pt gmon.out *.core *.log

.PHONY: all check clean
