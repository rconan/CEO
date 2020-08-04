include ../common.mk

all: $(texsrc)
	pdflatex ceo_manual.tex	
	makeindex ceo_manual
#	bibtex ceo_manual
	pdflatex ceo_manual.tex	
	touch ceo_manual.pdf
