#!/bin/sh 
pdflatex -shell-escape -synctex=1 -interaction=nonstopmode GA.tex
pdflatex -shell-escape -synctex=1 -interaction=nonstopmode GB.tex
lualatex -shell-escape  -interaction=nonstopmode dcs-lp.tex
clean_latex 	
