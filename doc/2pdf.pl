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


`dvips -P pdf -o tao_manual.ps tao_manual`; 
`ps2pdf tao_manual.ps`;
`mv tao_manual.pdf tao_manual-$rev.pdf`; 
`dvips -o tao_manual.ps tao_manual`;
`mv tao_manual.ps tao_manual-$rev.ps`; 
`cp tao_manual-*.p*  /home/dcs/public_html/bmad`;

open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");
$file = "/home/dcs/public_html/bmad/tao.html";
open (FM, $file) || die ("Cannot open File: $file\n");
while (<FM>) {
  if (/tao_manual-(\S+?)\.p+/) {s/$1/$rev/g;}
  print (F_OUT);
}

close (FM);
close (F_OUT);
`mv temp.out $file`;
