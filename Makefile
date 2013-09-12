include common.mk

SOURCE_DIR	= ceo source atmosphere imaging centroiding aaStats paStats iterativeSolvers

all: makefile 
	mkdir -p include lib
	for i in $(SOURCE_DIR); do (make -C $$i all); done

tex: $(texsrc)
	for i in $(SOURCE_DIR); do (make -C $$i tex); done
	rm -f doc/ceo.manual.main.tex
	for i in $(SOURCE_DIR); do (echo -e "\include{ceo.manual.$$i}\n">>doc/ceo.manual.main.tex); done
	for i in $(SOURCE_DIR); do (echo -e "\section{$$i}\n\label{sec:$$i}\n\n\input{../$$i/$$i}">doc/ceo.manual.$$i.tex); done

doc: tex
	make -C doc all

touch: 
	find . -name \*.nw -exec touch {} \;

makefile:
	for i in $(SOURCE_DIR); do (cp Makefile.common $$i/Makefile; sed -i -e "s/filename/$$i/g" $$i/Makefile); done

clean_makefile:
	for i in $(SOURCE_DIR); do (rm -f $$i/Makefile); done

clean:
	for i in $(SOURCE_DIR); do (make -C $$i clean); done
	rm -f *.*~

