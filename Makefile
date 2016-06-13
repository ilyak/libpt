CC= clang
CFLAGS= -Weverything -Wno-conversion -Wno-format-nonliteral -Wno-unused-parameter -fcolor-diagnostics -g
LDFLAGS=
LIBS= -lm

ALL_O= log.o vis.o xmalloc.o
RM= rm -f

all: pt

pt: pt.o $(ALL_O)
	$(CC) -o $@ $(CFLAGS) pt.o $(ALL_O) $(LDFLAGS) $(LIBS)

check: pt
	@./pt && echo success

clean:
	$(RM) $(ALL_O) pt pt.o *.core *.log

.PHONY: all check clean
