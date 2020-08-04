include ../common.mk

all: ActiveOptics.nw
	$(TANGLE) -RLinearActiveOptics.py  ActiveOptics.nw > LinearActiveOptics.py

doc: ActiveOptics.nw
	$(WEAVE) -delay -index  ActiveOptics.nw > LinearActiveOptics.tex
	pdflatex main.tex > log
	pdflatex main.tex > log
