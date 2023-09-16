#!/usr/bin/perl

`pdflatex tao`;
`makeindex tao`;

open (FC, "tao.ind") || die ("Cannot open File: tao.idx\n");
open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");

while (<FC>) {
  s/_/\\_/g;
  s/\%/\\\%/g;
  print (F_OUT);
}

close (FC);
close (F_OUT);
`mv temp.out tao.ind`;

`pdflatex tao`;
