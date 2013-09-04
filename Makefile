include common.mk

SOURCE_DIR	= ceo source atmosphere imaging centroiding

all: makefile
	mkdir -p include lib
	for i in $(SOURCE_DIR); do (make -C $$i all); done

tex: $(texsrc)
	for i in $(SOURCE_DIR); do (make -C $$i tex); done

doc: tex
	make -C doc all

touch: 
	find . -name \*.nw -exec touch {} \;

makefile:
	for i in $(SOURCE_DIR); do (cp Makefile.common $$i/Makefile; sed -i -e "s/filename/$$i/g" $$i/Makefile); done

clean:
	for i in $(SOURCE_DIR); do (make -C $$i clean); done
	rm -f *.*~

