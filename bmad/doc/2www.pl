#!/usr/bin/perl

# Make index

`./makeindex.pl`;

# To web page

$found = 0;
open (FC, "cover-page.tex") || die ("Cannot open File: cover-page.tex\n");
while (<FC>) {
  if (/Revision: +(.+) +\\\\/) {
    $date = $1;
    print "Revision: $date\n";
    $found = 1;
    last;
  }
}

close (FC);
if (! $found) {die ("Revision line not found in: cover_page.tex\n");}

`pdflatex bmad`; 
`chmod g+w bmad.pdf`;

open (F_IN, "manual_template.html") || die ("Cannot open File: manual_template.html\n");

open (F_OUT, ">basic_manual.html") || die ("Cannot open basic_manual.html file\n");

%month2num = qw(
    January 01  February 02  March 03  April 04  May 05  June 06
    July 07  August 08  September 09  October 10 November 11 December 12
);

$date =~ /(.+) (.+), (.+) */;
$day = $2;
if ($day < 10) { $day = "0$day"; }    # Two digit day

$rev = "$3-$month2num{$1}-$day";
print "Manual: bmad-manual-$rev\n";

while (<F_IN>) {
  s/RRR/$rev/g;
  s/DDD/$date/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`mv basic_manual.html  /nfs/classe/www/html/bmad/basic_manual.html`;

`cp bmad.pdf /nfs/classe/www/html/bmad/bmad-manual-$rev.pdf`;
`chmod g+w  /nfs/classe/www/html/bmad/bmad-manual-$rev.pdf`;

`cp bmad.pdf /nfs/classe/www/html/bmad/bmad-manual.pdf`;   # For Chris
