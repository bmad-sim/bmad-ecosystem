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

`pdflatex bmad`; 
`cp bmad.pdf bmad-manual-$rev.pdf`; 

`cp /home/dcs/public_html/bmad/manual.html .`;
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

`cp /home/dcs/public_html/bmad/manual.html /home/dcs/public_html/bmad/manual.html.$rev`;
`mv manual.html  /home/dcs/public_html/bmad`;

`cp bmad-manual-*.pdf ~/public_html/bmad`;
`cp bmad.pdf          ~/public_html/bmad/bmad-manual.pdf`;   # For Chris
