#!/usr/bin/perl

`latex bmad-manual`;
`makeindex bmad-manual`;

open (FC, "bmad-manual.ind") || die ("Cannot open File: bmad-manual.idx\n");
open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");

while (<FC>) {
  s/_/\\_/g;
  s/\%/\\\%/g;
  print (F_OUT);
}

close (FC);
close (F_OUT);
`mv temp.out bmad-manual.ind`;


`latex bmad-manual`;
