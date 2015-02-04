CEOPATH	        = /home/ubuntu/Dropbox/CEO
CUDAPATH	= /usr/local/cuda
PYTHONPATH      = /home/ubuntu/anaconda
CEOPYPATH	= $(CEOPATH)/python/ceo
NVCC          	= $(CUDAPATH)/bin/nvcc
CUDALIBPATH   	= $(CUDAPATH)/lib64
MATLABINCS	= -I/priv/monarcas1/rconan/MATLAB/R2013a/extern/include \
	-I/export/monarcas1/rconan/MATLAB/R2013a/toolbox/distcomp/gpu/extern/include
CUDALIBS	= cusparse cufft cublas cudart cuda 
NOWEBPATH	= /usr
WEAVE   	= $(NOWEBPATH)/bin/noweave
TANGLE    	= $(NOWEBPATH)/bin/notangle
CPIF	    	= $(NOWEBPATH)/bin/cpif
TEXTOPDF  	= pdflatex
NVCCFLAGS	= -gencode=arch=compute_20,code=\"sm_20,compute_20\" \
		--compiler-options=-ansi,-D_GNU_SOURCE,-fwrapv,-fPIC,-fno-omit-frame-pointer,-pthread,-fno-strict-aliasing -O2
LIBS 		= -L$(CEOPATH)/lib $(CUDALIBPATH:%=-L%) -lceo -lcurl -ljsmn $(CUDALIBS:%=-l%)
INCS		= -I. -I$(CEOPATH)/include #$(MATLABINCS)
SHELL		= /bin/bash

texsrc = $(nwsrc:%.nw=%.tex)
header = $(nwsrc:%.nw=%.h)
obj    = $(nwsrc:%.nw=%.o)
cusrc  = $(nwsrc:.%nw=%.cu)
libsrc = $(CEOPATH)/lib/libceo.a 

.SUFFIXES: .nw .tex .cu .mex .bin .py

.cu.o: 
	$(NVCC) $(INCS) $(NVCCFLAGS) -o $@ -c $<

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

.nw.mex:
	$(TANGLE) -L -R$@ $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	mv $@ $@.cu

.nw.bin:
	$(TANGLE) -L -R$@ $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	mv $@ $@.cu
	make -C $(CEOPATH) all
	$(NVCC) $(INCS) $(LIBS) $@.cu

.nw.py:
	$(TANGLE) -R$@ $< > $@
