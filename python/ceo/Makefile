include ../../common.mk

.pyx.so:
	cython --directive embedsignature=True --cplus $< -o  $<.cu
	$(NVCC) -w $(INCS) -I/home/ubuntu/anaconda/envs/python34/include/python3.4m -I/home/ubuntu/anaconda/envs/python34/lib/python3.4/site-packages/numpy/core/include -arch=sm_20 -lineinfo -Xcompiler '-fPIC,-pthread,-O3' -o  $<.o -c  $<.cu
	$(NVCC) -shared $<.o -o $@ $(LIBS)

aco:
	$(TANGLE) -RLinearActiveOptics.py  ActiveOptics.nw > LinearActiveOptics.py
	$(WEAVE) -delay -index  ActiveOptics.nw > doc/LinearActiveOptics.tex
	make -C doc
