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

`pdflatex tao`;

$file = "tao_template.html";
open (F_IN, $file) || die ("Cannot open File: $file\n");

open (F_OUT, ">basic_tao.html") || die ("Cannot open basic_tao.html file\n");

while (<F_IN>) {
  s/RRR/$rev/g;
  s/DDD/$date/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`mv basic_tao.html  /home/dcs/public_html/bmad`;
`cp tao.pdf   /home/dcs/public_html/bmad/tao-manual-$rev.pdf`;

