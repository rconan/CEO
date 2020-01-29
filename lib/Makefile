include ../common.mk

all:
	rm -f link.o
	$(NVCC) --device-link *.o -o link.o
	$(NVCC) --lib *.o -o libceo.a
