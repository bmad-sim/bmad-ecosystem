#!/usr/bin/perl

`latex tao_manual`;
`makeindex tao_manual`;

open (FC, "tao_manual.ind") || die ("Cannot open File: tao_manual.idx\n");
open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");

while (<FC>) {
  s/_/\\_/g;
  s/\%/\\\%/g;
  print (F_OUT);
}

close (FC);
close (F_OUT);
`mv temp.out tao_manual.ind`;


`latex tao_manual`;
