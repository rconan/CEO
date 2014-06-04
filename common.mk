CEOPATH	        = /home/rconan/CEO
CUDAPATH	= /opt/local/cuda
NVCC          	= $(CUDAPATH)/bin/nvcc
CUDALIBPATH   	= $(CUDAPATH)/lib64
MATLABINCS	= -I/priv/monarcas1/rconan/MATLAB/R2013a/extern/include \
	-I/export/monarcas1/rconan/MATLAB/R2013a/toolbox/distcomp/gpu/extern/include
CUDALIBS	= cusparse cufft cublas cudart cuda 
NOWEBPATH	= /opt/local/noweb
WEAVE   	= $(NOWEBPATH)/bin/noweave
TANGLE    	= $(NOWEBPATH)/bin/notangle
CPIF	    	= $(NOWEBPATH)/bin/cpif
TEXTOPDF  	= pdflatex
NVCCFLAGS	= -gencode=arch=compute_20,code=\"sm_20,compute_20\" \
		--compiler-options=-ansi,-D_GNU_SOURCE,-fPIC,-fno-omit-frame-pointer,-pthread -O2
LIBS 		= -L$(CEOPATH)/lib $(CUDALIBPATH:%=-L%) $(CUDALIBS:%=-l%)
INCS		= -I. -I$(CEOPATH)/include #$(MATLABINCS)

texsrc = $(nwsrc:%.nw=%.tex)
header = $(nwsrc:%.nw=%.h)
obj    = $(nwsrc:%.nw=%.o)
cusrc  = $(nwsrc:.%nw=%.cu)
libsrc = lib/libceo.a 

.SUFFIXES: .nw .tex .cu .mex .bin

.cu.o:
	$(NVCC) $(INCS) $(NVCCFLAGS) -o $@ -c $<

.nw.tex:
	$(WEAVE) -delay -index $< > $@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	sed -i -e "s/label{eq:/label{$*.eq:/g" -e "s/ref{eq:/ref{$*.eq:/g" $@
	sed -i -e "s/label{fig:/label{$*.fig:/g" -e "s/ref{fig:/ref{$*.fig:/g" $@
	sed -i -e "s/label{sec:/label{$*.sec:/g" -e "s/ref{sec:/ref{$*.sec:/g" $@
.nw.h:
	$(TANGLE) -R$@ $< > $@  
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $@
	mv $@ ../include/
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
	$(NVCC) $(INCS) $(LIBS) $@.cu -lceo -lcurl -ljsmn
