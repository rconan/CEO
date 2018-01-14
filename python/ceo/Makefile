include ../../common.mk

.pyx.so:
	cython --directive embedsignature=True --cplus $< -o  $<.cu
	$(NVCC) -w $(INCS) -I/home/ubuntu/anaconda/include/python3.6m -I/home/ubuntu/anaconda/lib/python3.6/site-packages/numpy/core/include -lineinfo -Xcompiler '-fPIC,-pthread,-O3' -o  $<.o -c  $<.cu
	$(NVCC) -shared $<.o -o $@ $(LIBS)

aco:
	$(TANGLE) -RLinearActiveOptics.py  ActiveOptics.nw > LinearActiveOptics.py
	$(WEAVE) -delay -index  ActiveOptics.nw > doc/LinearActiveOptics.tex
	make -C doc
