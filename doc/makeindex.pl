#!/usr/bin/perl

`latex bmad`;
`makeindex bmad`;
`makeindex bmad.rdx -o bmad.rnd`;

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
