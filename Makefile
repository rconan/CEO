CUDAPATH	= /opt/local/cuda
NVCC          	= $(CUDAPATH)/bin/nvcc
CUDALIBPATH   	= $(CUDAPATH)/lib64
CUDALIBS	= cufft cudart cuda
NOWEBPATH	= /opt/local/noweb
WEAVE   	= $(NOWEBPATH)/bin/noweave
TANGLE    	= $(NOWEBPATH)/bin/notangle
CPIF	    	= $(NOWEBPATH)/bin/cpif
TEXTOPDF  	= pdflatex
NVCCFLAGS	= -gencode=arch=compute_20,code=\"sm_20,compute_20\" \
		--compiler-options=-ansi,-D_GNU_SOURCE,-fPIC,-fno-omit-frame-pointer,-pthread -O2
LIBS 		= -L./lib $(CUDALIBPATH:%=-L%) $(CUDALIBS:%=-l%)
CDIR		= ./imaging
INCS		= -I$(CDIR) -I./include -I/priv/monarcas1/rconan/MATLAB/R2013a/extern/include \
		 -I/export/monarcas1/rconan/MATLAB/R2013a/toolbox/distcomp/gpu/extern/include

nwsrc = ceo.nw source.nw atmosphere.nw imaging.nw centroiding.nw
texsrc = $(nwsrc:%.nw=%.tex)
header = $(nwsrc:%.nw=%.h)
obj    = $(nwsrc:%.nw=%.o)
cusrc  = $(nwsrc:.%nw=%.cu)

all: lib tex

lib: $(nwsrc) $(header) $(cusrc) $(obj)
	ar rcs libceo.a $(obj)
	mv libceo.a ./lib

tex: $(texsrc)

clean:
	rm -f $(obj)

.SUFFIXES: .nw .tex .cu .mex

.cu.o: $(header)
	$(NVCC) $(INCS) $(NVCCFLAGS) -o $@ -c ./src/$<

.nw.tex:
	$(WEAVE) -index $< > ./doc/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ./doc/$@
	$(TEXTOPDF) -output-directory doc ./doc/$@
	$(TEXTOPDF) -output-directory doc ./doc/$@
.nw.h:
	$(TANGLE) -R$@ $< > ./include/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ./include/$@
.nw.cu:
	$(TANGLE) -L -R$@ $< > ./src/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ./src/$@

.nw.mex:
	$(TANGLE) -R$@ $< > $(CDIR)/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $(CDIR)/$@
	mv $(CDIR)/$@ $(CDIR)/$@.cu

