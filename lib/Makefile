include ../common.mk

all: libceo.a *.o

libceo.a: 
#	touch link.o
	$(NVCC) $(NVCCFLAGS)  --device-link *.o -o link.o -lcudadevrt -lcudart
	$(NVCC) --lib *.o -o libceo.a
