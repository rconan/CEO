include ../../common.mk

.pyx.so:
	cython --directive embedsignature=True --directive language_level=2 --cplus $< -o  $<.cu
	$(NVCC) -w $(INCS) -I$(PYTHONPATH)/include/python3.7m -I$(PYTHONPATH)/lib/python3.7/site-packages/numpy/core/include -lineinfo -Xcompiler '-fPIC,-pthread,-O3' -o  $<.o -c  $<.cu
	$(NVCC) -shared $<.o -o $@ $(LIBS)

aco:
	$(TANGLE) -RLinearActiveOptics.py  ActiveOptics.nw > LinearActiveOptics.py
	$(WEAVE) -delay -index  ActiveOptics.nw > doc/LinearActiveOptics.tex
	make -C doc
