#!/usr/bin/perl

$found = 0;
open (FC, "cover-page.tex") || die ("Cannot open File: cover-page.tex\n");
while (<FC>) {
  if (/Revision: +(\S+) +\\\\/) {
    $rev = $1;
    print "Revision: $rev\n";
    $found = 1;
  }
}

close (FC);
if (! $found) {die ("Revision line not found in: cover_page.tex\n");}


`latex bmad-manual`; 
`dvips -o bmad-manual.ps -P generic bmad-manual`; 
`ps2pdf bmad-manual.ps`;
`rm /home/dcs/public_html/bmad/bmad-manual-*pdf`;
`rm /home/dcs/public_html/bmad/bmad-manual-*ps`;
`mv bmad-manual.pdf /home/dcs/public_html/bmad/bmad-manual-$rev.pdf`; 
`dvips -o bmad-manual.ps bmad-manual`;
`mv bmad-manual.ps /home/dcs/public_html/bmad/bmad-manual-$rev.ps`; 

open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");
$file = "/home/dcs/public_html/bmad/manual.html";
open (FM, $file) || die ("Cannot open File: $file\n");
while (<FM>) {
  if (/bmad-manual-(\S+?)\.p+/) {s/$1/$rev/g;}
  print (F_OUT);
}

close (FM);
close (F_OUT);
`mv temp.out $file`;
