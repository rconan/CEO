include ../../common.mk

PYV=$(shell python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)")

.pyx.so:
	cython -2 --directive embedsignature=True --cplus $< -o  $<.cu
	$(NVCC) -w $(INCS) -I$(PYTHONPATH)/include/python$(PYV) -I$(PYTHONPATH)/lib/python$(PYV)/site-packages/numpy/core/include -lineinfo -Xcompiler '-fPIC,-pthread,-O3' -o  $<.o -c  $<.cu
	$(NVCC) -shared $<.o -o $@ $(LIBS)

aco:
	$(TANGLE) -RLinearActiveOptics.py  ActiveOptics.nw > LinearActiveOptics.py
	$(WEAVE) -delay -index  ActiveOptics.nw > doc/LinearActiveOptics.tex
	make -C doc
