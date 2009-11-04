#!/usr/bin/perl

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

`latex tao`;
`dvips -P pdf -o tao.ps -z tao`; 
`ps2pdf tao.ps`;
`rm tao-manual-*pdf`;
`rm tao-manual-*ps`;
`cp tao.pdf tao-manual-$rev.pdf`; 
`dvips -o tao.ps tao`;
`mv tao.ps tao-manual-$rev.ps`; 

open (F_OUT, ">tao.html") || die ("Cannot open tao.html file\n");
$file = "tao_template.html";
open (F_IN, $file) || die ("Cannot open File: $file\n");
while (<F_IN>) {
  s/RRR/$rev/g;
  s/DDD/$date/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`scp tao.html       dcs\@lnx209.lepp.cornell.edu:/home/dcs/public_html/bmad`;
`scp tao-manual-*.pdf dcs\@lnx209.lepp.cornell.edu:/home/dcs/public_html/bmad`;
`scp tao-manual-*.ps  dcs\@lnx209.lepp.cornell.edu:/home/dcs/public_html/bmad`;

