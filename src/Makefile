LIBPT_A= libpt.a
LIBPT_O= pt.o

AR= ar rc
RANLIB= ranlib

$(LIBPT_A): $(LIBPT_O)
	$(AR) $@ $(LIBPT_O)
	$(RANLIB) $@

clean:
	rm -f $(LIBPT_A) $(LIBPT_O)

.PHONY: clean
