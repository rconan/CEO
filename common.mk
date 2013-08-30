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
