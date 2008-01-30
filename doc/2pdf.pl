#!/usr/bin/perl

`scp lnx209:/home/dcs/public_html/bmad/tao.html .`;

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

`latex tao`;
`dvips -P pdf -o tao.ps -z tao`; 
`ps2pdf tao.ps`;
`mv tao.pdf tao-manual-$rev.pdf`; 
`dvips -o tao.ps tao`;
`mv tao.ps tao-manual-$rev.ps`; 

open (F_OUT, ">temp.out") || die ("Cannot open temparary file\n");
$file = "tao.html";
open (FM, $file) || die ("Cannot open File: $file\n");
while (<FM>) {
  if (/tao-manual-(\S+?)\.p+/) {s/$1/$rev/g;}
  print (F_OUT);
}

close (FM);
close (F_OUT);
`mv temp.out $file`;

`scp tao.html       lnx209:/home/dcs/public_html/bmad`;
`scp tao-manual-*.pdf lnx209:/home/dcs/public_html/bmad`;
`scp tao-manual-*.ps  lnx209:/home/dcs/public_html/bmad`;

