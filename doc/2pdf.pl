#!/usr/bin/perl

`scp dcs\@lnx209.lns.cornell.edu:/home/dcs/public_html/bmad/manual.html .`;

$found = 0;
open (FC, "cover-page.tex") || die ("Cannot open File: cover-page.tex\n");
while (<FC>) {
  if (/Revision: +(\S+) +\\\\/) {
    $rev = $1;
    print "Revision: $rev\n";
    $found = 1;
    $date_line = <FC>;
    $date_line =~ /^ *(\S.+\S) +\\\\/;
    $date = $1;
    print "Date: $date\n";
    last;
  }
}

close (FC);
if (! $found) {die ("Revision line not found in: cover_page.tex\n");}

`latex bmad`; 
`dvips -o bmad.ps -z bmad`; 
`ps2pdf bmad.ps`;
`rm bmad-manual-*pdf`;
`rm bmad-manual-*ps`;
`mv bmad.pdf bmad-manual-$rev.pdf`; 
`dvips -o bmad.ps bmad`;
`mv bmad.ps bmad-manual-$rev.ps`; 

open (F_OUT, ">manual.html") || die ("Cannot open manual.html file\n");
$file = "manual_template.html";
open (F_IN, $file) || die ("Cannot open File: $file\n");
while (<F_IN>) {
  s/RRR/$rev/g;
  s/DDD/$date/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`scp manual.html       dcs\@lnx209.lns.cornell.edu:/home/dcs/public_html/bmad`;
`scp bmad-manual-*.pdf dcs\@lnx209.lns.cornell.edu:/home/dcs/public_html/bmad`;
`scp bmad-manual-*.ps  dcs\@lnx209.lns.cornell.edu:/home/dcs/public_html/bmad`;
