#!/bin/bash

# Do this on a system with pdfcrop 
for f in *ps
  do
    base="${f%.*}"
    echo $f '-> ' $base.pdf
    #sed  's/5 mul /1 mul /' -i '' $f
    sed  's/5 mul /1 mul /' -i $f
    ps2pdf $f
    pdfcrop $base.pdf $base-crop.pdf
    mv $base-crop.pdf $base.pdf
  done
