#!/usr/bin/perl

`latex bmad`;
`makeindex bmad`;

#open (FC, "bmad.ind") || die ("Cannot open File: bmad.idx\n");
#open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");
#while (<FC>) {
#  s/_/\\_/g;
#  s/\%/\\\%/g;
#  print (F_OUT);
#}
#close (FC);
#close (F_OUT);
#`mv temp.out bmad.ind`;

`latex bmad`;
