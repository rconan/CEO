include ../common.mk

nwsrc = centroiding.nw

all: lib tex

lib:  $(header) $(obj)
	$(AR) rcs ../$(libsrc) $(obj)

tex: $(texsrc)

clean:
	rm -f $(obj) *.*~

