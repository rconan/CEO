include common.mk

nwsrc = ceo.nw source.nw atmosphere.nw imaging.nw centroiding.nw
texsrc = $(nwsrc:%.nw=%.tex)
header = $(nwsrc:%.nw=%.h)
obj    = $(nwsrc:%.nw=%.o)
cusrc  = $(nwsrc:.%nw=%.cu)

all: lib tex

lib: $(nwsrc) $(header) $(cusrc) $(obj)
	ar rcs libceo.a $(obj)
	mv libceo.a ./lib

tex: $(texsrc)

clean:
	rm -f $(obj)

.SUFFIXES: .nw .tex .cu .mex

.cu.o: $(header)
	$(NVCC) $(INCS) $(NVCCFLAGS) -o $@ -c ./src/$<

.nw.tex:
	$(WEAVE) -index $< > ./doc/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ./doc/$@
	$(TEXTOPDF) -output-directory doc ./doc/$@
	$(TEXTOPDF) -output-directory doc ./doc/$@
.nw.h:
	$(TANGLE) -R$@ $< > ./include/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ./include/$@
.nw.cu:
	$(TANGLE) -L -R$@ $< > ./src/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' ./src/$@

.nw.mex:
	$(TANGLE) -L -R$@ $< > $(CDIR)/$@
	sed -i -e 's/LLL/<<</g' -e 's/RRR/>>>/g' $(CDIR)/$@
	mv $(CDIR)/$@ $(CDIR)/$@.cu

