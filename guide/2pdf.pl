#!/usr/bin/perl

$found = 0;
open (FC, "cover_page.tex") || die ("Cannot open File: cover_page.tex\n");
while (<FC>) {
  if (/Revision: +(\S+) +\\\\/) {
    $rev = $1;
    print "Revision: $rev\n";
    $found = 1;
  }
}

close (FC);
if (! $found) {die ("Revision line not found in: cover_page.tex\n");}


`dvips -P pdf bmad_manual`; 
`ps2pdf bmad_manual.ps`;
`mv bmad_manual.pdf bmad_manual-$rev.pdf`; 
`dvips bmad_manual`;
`mv bmad_manual.ps bmad_manual-$rev.ps`; 
`cp bmad_manual-*.p*  /home/dcs/public_html/bmad`;

open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");
$file = "/home/dcs/public_html/bmad/manual.html";
open (FM, $file) || die ("Cannot open File: $file\n");
while (<FM>) {
  if (/bmad_manual-(\S+?)\.p+/) {s/$1/$rev/g;}
  print (F_OUT);
}

close (FM);
close (F_OUT);
`mv temp.out $file`;
