include ../common.mk

nwsrc = source.nw

all: lib tex

lib:  $(header) $(obj)
	$(AR) rcs ../$(libsrc) $(obj)

tex: $(texsrc)

clean:
	rm -f $(obj) *.*~

