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
`chmod g+w bmad.pdf`;
`mv bmad.pdf bmad-manual-$rev.pdf`; 

open (F_IN, "manual_template.html") || die ("Cannot open File: manual_template.html\n");

open (F_OUT, ">manual.html") || die ("Cannot open manual.html file\n");

while (<F_IN>) {
  s/RRR/$rev/g;
  s/DDD/$date/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`cp manual.html  ~/public_html/bmad/manual.html`;

`cp bmad-manual-$rev.pdf ~/public_html/bmad`;
`chmod g+w  ~/public_html/bmad/bmad-manual-$rev.pdf`;

`cp bmad-manual-$rev.pdf ~/public_html/bmad/bmad-manual.pdf`;   # For Chris
