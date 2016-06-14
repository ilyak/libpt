CC= clang
CFLAGS= -I../libxutil -Weverything -Wno-conversion -Wno-format-nonliteral -Wno-unused-parameter -fcolor-diagnostics -g
LDFLAGS= -L../libxutil
LIBS= -lm -lxutil

ALL_O=
RM= rm -f

all: pt

pt: pt.o $(ALL_O)
	$(CC) -o $@ $(CFLAGS) pt.o $(ALL_O) $(LDFLAGS) $(LIBS)

check: pt
	@./pt -t h2o.dat && echo success
	@./pt -t nh3.dat && echo success
	@./pt -t ch4.dat && echo success

clean:
	$(RM) $(ALL_O) pt pt.o *.core *.log

.PHONY: all check clean
