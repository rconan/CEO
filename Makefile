include common.mk

#ls -d */ | sed -e 's,//$,,' -e 's,doc,,' -e 's,lib,,'  -e 's,include,,' | xargs
SOURCE_DIR	= utilities source atmosphere imaging centroiding shackHartmann aaStats BTBT GBTBT iterativeSolvers LMMSE plotly
TUTORIAL	= lgsao ngsao ltao ltaoVsAst geaos
PYTHON_DIR	= source atmosphere

all: makefile 
	mkdir -p include lib
	for i in $(SOURCE_DIR); do (make -C $$i all);echo -e "\n"; done

tex: $(texsrc)
	for i in $(SOURCE_DIR); do (make -C $$i tex); done
	for i in $(TUTORIAL); do (make -C TUTORIAL $$i.tex); done
	rm -f doc/ceo.manual.main.tex
	for i in $(SOURCE_DIR); do (echo -e "\input{ceo.manual.$$i}\n">>doc/ceo.manual.main.tex); done
	for i in $(SOURCE_DIR); do (echo -e "\section{$$i}\n\label{sec:$$i}\n\n\input{../$$i/$$i}">doc/ceo.manual.$$i.tex); done

cython:
	mkdir -p $(CEOPATH)/python
	rm -f $(CEOPATH)/python/ceo.pxd $(CEOPATH)/python/ceo.pyx
	cp $(CEOPATH)/etc/ceo.pxd $(CEOPATH)/python/ceo.pxd
	cp $(CEOPATH)/etc/ceo.pyx $(CEOPATH)/python/ceo.pyx
	for i in $(PYTHON_DIR); do (make -C $$i python);echo -e "\n"; done
	cython --cplus $(CEOPATH)/python/ceo.pyx -o $(CEOPATH)/python/ceo.cu
	$(NVCC) $(INCS) -I/home/rconan/anaconda/include/python2.7/ -I/home/rconan/anaconda/pkgs/numpy-1.8.1-py27_0/lib/python2.7/site-packages/numpy/core/include $(NVCCFLAGS) -o $(CEOPATH)/python/ceo.o -c $(CEOPATH)/python/ceo.cu
	$(NVCC) $(LIBS) -shared $(CEOPATH)/python/ceo.o -o $(CEOPATH)/python/ceo.so -lceo -lcurl -ljsmn
	rm -f $(CEOPATH)/python/ceo.cu $(CEOPATH)/python/ceo.o

doc: tex
	make -C doc all

touch: 
	find . -name \*.nw -exec touch {} \;

makefile: Makefile.common
	for i in $(SOURCE_DIR); do (cp Makefile.common $$i/Makefile; sed -i -e "s/filename/$$i/g" $$i/Makefile); done

jsmnlib: 
	make -C jsmn
	cp -P jsmn/jsmn.h include/
	cp -P jsmn/libjsmn.a lib/


clean_makefile:
	for i in $(SOURCE_DIR); do (rm -f $$i/Makefile); done

clean:
	for i in $(SOURCE_DIR); do (make -C $$i clean); done
	rm -f *.*~
	rm -f lib/libceo.a

cleanbins: makefile
	for i in $(SOURCE_DIR); do (make -C $$i cleanbins); done
