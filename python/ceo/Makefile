include ../../common.mk

.pyx.so:
	cython --cplus $< -o  $<.cu
	$(NVCC) $(INCS) -I$(PYTHONPATH)/include/python2.7/ -I$(PYTHONPATH)/lib/python2.7/site-packages/numpy/core/include $(NVCCFLAGS) -o  $<.o -c  $<.cu
	$(NVCC) -shared $<.o -o $@ $(LIBS) 
