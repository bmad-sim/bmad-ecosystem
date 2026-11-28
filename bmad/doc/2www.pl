#!/usr/bin/perl

# Make index

`./makeindex.pl`;
`git pull`;

# To web page

$thisdate = '';
open (FC, "cover-page.tex") || die ("Cannot open File: cover-page.tex\n");
while (<FC>) {
  if (/Date: +(\S+) /) {
    $thisdate = $1;
    print "Date: $thisdate\n";
    last;
  }
}

close (FC);
if (! $thisdate) {die ("Date not found in: cover-page.tex\n");}

`pdflatex bmad`; 
`chmod g+w bmad.pdf`;

open (F_IN, "manual_template.html") || die ("Cannot open File: manual_template.html\n");
open (F_OUT, ">basic_manual.html") || die ("Cannot open basic_manual.html file\n");

print "Manual: bmad-manual-$thisdate\n";

while (<F_IN>) {
  s/DDD/$thisdate/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`mv basic_manual.html  /nfs/classe/www/html/bmad/basic_manual.html`;

`cp bmad.pdf /nfs/classe/www/html/bmad/bmad-manual-$thisdate.pdf`;
`chmod g+w  /nfs/classe/www/html/bmad/bmad-manual-$thisdate.pdf`;

`cp bmad.pdf /nfs/classe/www/html/bmad/bmad-manual.pdf`;   # For Chris
