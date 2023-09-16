#!/bin/bash

f2py -c geodesiclm.pyf -L../geodesicLM/ -lgeodesiclm --fcompiler=gnu95 -L/usr/local/lib -llapack -L/usr/local/lib -lblas -lgfortran
