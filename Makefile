include common.mk

#ls -d */ | sed -e 's,//$,,' -e 's,doc,,' -e 's,lib,,'  -e 's,include,,' | xargs
SOURCE_DIR	= utilities source atmosphere imaging centroiding shackHartmann aaStats BTBT GBTBT iterativeSolvers LMMSE plotly system
TUTORIAL	= ngsao lgsao ltao ltaoVsAst geaos
PYTHON_DIR	= utilities source atmosphere centroiding imaging shackHartmann aaStats GBTBT LMMSE

all: makefile jsmnlib
	mkdir -p include lib
	for i in $(SOURCE_DIR); do (make -C $$i all);echo -e "\n"; done

tex: makefile $(texsrc)
	for i in $(SOURCE_DIR); do (make -C $$i tex); done
	for i in $(TUTORIAL); do (make -C TUTORIAL $$i.tex); done
	rm -f doc/ceo.manual.main.tex
	for i in $(SOURCE_DIR); do (echo -e "\input{ceo.manual.$$i}\n">>doc/ceo.manual.main.tex); done
	for i in $(SOURCE_DIR); do (echo -n "\chapter" >doc/ceo.manual.$$i.tex; echo -e "{$$i}\n\label{sec:$$i}\n\n\input{../$$i/$$i}">>doc/ceo.manual.$$i.tex); done

cython: all 
	mkdir -p $(CEOPYPATH)
	rm -f $(CEOPYPATH)/ceo.pxd $(CEOPYPATH)/ceo.pyx
	cp $(CEOPATH)/etc/ceo.pxd $(CEOPYPATH)/ceo.pxd
	cp $(CEOPATH)/etc/ceo.pyx $(CEOPYPATH)/ceo.pyx
	for i in $(PYTHON_DIR); do (echo -e "\n# $$i.nw">>$(CEOPYPATH)/ceo.pxd;echo -e "\n# $$i.nw">>$(CEOPYPATH)/ceo.pyx;make -C $$i python);echo -e "\n"; done
	cython --cplus $(CEOPYPATH)/ceo.pyx -o $(CEOPYPATH)/ceo.cu
	$(NVCC) $(INCS) -I$(PYTHONPATH)/include/python2.7/ -I$(PYTHONPATH)/lib/python2.7/site-packages/numpy/core/include $(NVCCFLAGS) -o $(CEOPYPATH)/ceo.o -c $(CEOPYPATH)/ceo.cu
	$(NVCC) $(LIBS) -shared $(CEOPYPATH)/ceo.o -o $(CEOPYPATH)/ceo.so -lceo -lcurl -ljsmn
	rm -f $(CEOPYPATH)/ceo.cu $(CEOPYPATH)/ceo.o
	set PYTHONPATH=$(CEOPYPATH)

doc: tex
	make -C doc all

touch: 
	find . -name \*.nw -exec touch {} \;

makefile: Makefile.common
	for i in $(SOURCE_DIR); do (cp Makefile.common $$i/Makefile; sed -i -e "s/filename/$$i/g" $$i/Makefile); done

jsmnlib: 
	mkdir -p include lib
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
