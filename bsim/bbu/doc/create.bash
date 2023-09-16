#!/bin/bash
texi2pdf bbu_doc.tex
rm -f *.bbl *.log *.aux *.blg bbu_doc.dvi bbu_doc.ps
