include ../../common.mk

.pyx.so:
	cython --directive embedsignature=True --cplus $< -o  $<.cu
	$(NVCC) -w $(INCS) -I$(PYTHONPATH)/include/python2.7/ -I$(PYTHONPATH)/lib/python2.7/site-packages/numpy/core/include -arch=sm_20 -lineinfo -Xcompiler '-fPIC,-pthread,-O3' -o  $<.o -c  $<.cu
	$(NVCC) -shared $<.o -o $@ $(LIBS)

aco:
	$(TANGLE) -RLinearActiveOptics.py  ActiveOptics.nw > LinearActiveOptics.py
	$(WEAVE) -delay -index  ActiveOptics.nw > doc/LinearActiveOptics.tex
	make -C doc
