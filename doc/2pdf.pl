#!/usr/bin/perl

`scp lnx209:/home/dcs/public_html/bmad/manual.html .`;

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


`latex bmad`; 
`dvips -o bmad.ps bmad`; 
`ps2pdf bmad.ps`;
`rm bmad-*pdf`;
`rm bmad-*ps`;
`mv bmad.pdf bmad-manual-$rev.pdf`; 
`dvips -o bmad.ps bmad`;
`mv bmad.ps bmad-manual-$rev.ps`; 

open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");
$file = "manual.html";
open (FM, $file) || die ("Cannot open File: $file\n");
while (<FM>) {
  if (/bmad-manual-(\S+?)\.p+/) {s/$1/$rev/g;}
  print (F_OUT);
}

close (FM);
close (F_OUT);
`mv temp.out $file`;

`scp manual.html       lnx209:/home/dcs/public_html/bmad`;
`scp bmad-manual-*.pdf lnx209:/home/dcs/public_html/bmad`;
`scp bmad-manual-*.ps  lnx209:/home/dcs/public_html/bmad`;
