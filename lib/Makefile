include ../common.mk

all:
	touch link.o
	$(NVCC) $(NVCCFLAGS)  --device-link *.o -o link.o
	$(NVCC) --lib *.o -o libceo.a
