#!/usr/bin/perl

`latex bmad_manual`;
`makeindex bmad_manual`;

open (FC, "bmad_manual.ind") || die ("Cannot open File: bmad_manual.idx\n");
open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");

while (<FC>) {
  s/_/\\_/g;
  s/\%/\\\%/g;
  print (F_OUT);
}

close (FC);
close (F_OUT);
`mv temp.out bmad_manual.ind`;


`latex bmad_manual`;
