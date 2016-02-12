include common.mk

#ls -d */ | sed -e 's,//$,,' -e 's,doc,,' -e 's,lib,,'  -e 's,include,,' | xargs
SOURCE_DIR	= utilities source atmosphere imaging centroiding shackHartmann aaStats BTBT GBTBT iterativeSolvers LMMSE plotly rayTracing gmtMirrors segmentPistonSensor pyramid
TUTORIAL	= ngsao lgsao ltao ltaoVsAst geaos
CYTHON_DIR	= utilities rayTracing source imaging centroiding shackHartmann atmosphere LMMSE aaStats gmtMirrors segmentPistonSensor pyramid

all: makefile jsmnlib
ifeq ($(wildcard include/plotly.credentials), )
	echo "plotly.credentials doesn't exist!"
	cp include/plotly.credentials.sample include/plotly.credentials
else
	echo "plotly.credentials does exist!"
endif
	mkdir -p include lib
	for i in $(SOURCE_DIR); do (make -C $$i src);echo -e "\n"; done
	for i in $(SOURCE_DIR); do (make -C $$i lib);echo -e "\n"; done
	make -C lib all

tex: makefile $(texsrc)
	for i in $(SOURCE_DIR); do (make -C $$i tex); done
	for i in $(TUTORIAL); do (make -C TUTORIAL $$i.tex); done
	rm -f doc/ceo.manual.main.tex
	for i in $(SOURCE_DIR); do (echo -e "\include{ceo.manual.$$i}\n">>doc/ceo.manual.main.tex); done
	for i in $(SOURCE_DIR); do (echo -n "\chapter" >doc/ceo.manual.$$i.tex; echo -e "{$$i}\n\label{sec:$$i}\n\n\input{../$$i/$$i}">>doc/ceo.manual.$$i.tex); done


cython: makefile
	for i in $(CYTHON_DIR); do (make -C $$i cysrc);echo -e "\n"; done
	for i in $(CYTHON_DIR); do (make -C $$i cylib);echo -e "\n"; done

doc: tex
	make -C doc all

pydoc: cython
	make -C python/docs html 

test:
	make -C test all

ripython:
	env PYTHONPATH="$(CEOPATH)/python" ipython notebook --no-browser

ipythonserver:
	env PYTHONPATH="$(CEOPATH)/python" ipython notebook --profile=nbserver

ipython:
	env PYTHONPATH="$(CEOPATH)/python" ipython #notebook

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

cleanpython:
	for i in $(SOURCE_DIR); do (make -C $$i cleanpython); done
	rm -f python/ceo/*.so
	rm -f python/ceo/*.pxd
	rm -f python/ceo/*.pyx*

clean:
	for i in $(SOURCE_DIR); do (make -C $$i clean); done
	rm -f *.*~
	rm -f lib/libceo.a
	rm -f python/ceo/*.so
	rm -f python/ceo/*.pxd
	rm -f python/ceo/*.pyx*

cleanbins: makefile
	for i in $(SOURCE_DIR); do (make -C $$i cleanbins); done
