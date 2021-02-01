CEOPATH 	= $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
CUDAPATH	= /usr/local/cuda
PYTHONPATH      = $(HOME)/anaconda3
CEOPYPATH	= $(CEOPATH)/python/ceo
NVCC          	= $(CUDAPATH)/bin/nvcc
CUDALIBPATH   	= $(CUDAPATH)/lib64
CUDAINCPATH   	= $(CUDAPATH)/include
MATLABINCS	= -I/priv/monarcas1/rconan/MATLAB/R2013a/extern/include \
	-I/export/monarcas1/rconan/MATLAB/R2013a/toolbox/distcomp/gpu/extern/include
CUDALIBS	= curand cusparse cufft cublas cudart
NOWEBPATH	= /usr
WEAVE   	= $(NOWEBPATH)/bin/noweave
TANGLE    	= $(NOWEBPATH)/bin/notangle
CPIF	    	= $(NOWEBPATH)/bin/cpif
TEXTOPDF  	= pdflatex
NVCCFLAGS	= -DCEOPATH='"$(CEOPATH)"' -lineinfo -Xcompiler '-ansi,-D_GNU_SOURCE,-fwrapv,-fPIC,-fno-omit-frame-pointer,-pthread,-fno-strict-aliasing,-O3'
LIBS 		= -L$(CEOPATH)/lib $(CUDALIBPATH:%=-L%) -lceo -lcurl -ljsmn $(CUDALIBS:%=-l%)
INCS		= -I. -I$(CEOPATH)/include $(CUDAINCPATH:%=-I%) -I$(PYTHONPATH)/include #$(MATLABINCS)
SHELL		= /bin/bash

-include $(CEOPATH)/user.mk

texsrc = $(nwsrc:%.nw=%.tex)
header = $(nwsrc:%.nw=%.h)
obj    = $(nwsrc:%.nw=%.o)
sobj    = $(nwsrc:%.nw=%.so)
cusrc  = $(nwsrc:%.nw=%.cu)
nwcysrc  = $(nwsrc:%.nw=%.cython.nw)
texcysrc = $(nwcysrc:%.nw=%.tex)
pxdcysrc  = $(nwcysrc:%.nw=%.pxd)
pyxcysrc  = $(nwcysrc:%.nw=%.pyx)
pxdsrc  = $(nwcysrc:%.cython.nw=%.pxd)
pyxsrc  = $(nwcysrc:%.cython.nw=%.pyx)
libsrc = $(CEOPATH)/lib/libceo.a 

.SUFFIXES: .nw .tex .cu .mex .bin .py .pxd .pyx .so

.cu.o: 
	$(NVCC) $(INCS) $(NVCCFLAGS) --device-c -o $@ $<
	cp $@ ../lib/

.nw.tex:
	$(WEAVE) -delay -index $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	sed -i -e "s/label{eq:/label{$*.eq:/g" -e "s/ref{eq:/ref{$*.eq:/g" $@
	sed -i -e "s/label{fig:/label{$*.fig:/g" -e "s/ref{fig:/ref{$*.fig:/g" $@
	sed -i -e "s/label{tab:/label{$*.fig:/g" -e "s/ref{tab:/ref{$*.fig:/g" $@
	sed -i -e "s/label{sec:/label{$*.sec:/g" -e "s/ref{sec:/ref{$*.sec:/g" $@
.nw.h:
	$(TANGLE) -R$@ $< | $(CPIF) $@  
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	cp $@ ../include/
.nw.cu: 
	$(TANGLE) -L -R$@ $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	touch $(nwcysrc)

.nw.mex:
	$(TANGLE) -L -R$@ $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	mv $@ $@.cu

.nw.bin:
	$(TANGLE) -L -R$@ $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	mv $@ $@.cu
	make -C $(CEOPATH) all
	$(NVCC) $(NVCCFLAGS) -lineinfo $(INCS) $(LIBS) $@.cu

.nw.py:
	$(TANGLE) -R$@ $< > $@

.nw.pxd:
	$(TANGLE) -R$(pxdsrc) $< > $@
	cp $@ $(CEOPYPATH)/$(pxdsrc)

.nw.pyx:
	$(TANGLE) -R$(pyxsrc)  $< > $@
	cp $@ $(CEOPYPATH)/$(pyxsrc)

