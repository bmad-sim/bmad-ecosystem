#!/usr/bin/perl

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

`pdflatex tao`;

$file = "tao_template.html";
open (F_IN, $file) || die ("Cannot open File: $file\n");

open (F_OUT, ">basic_tao.html") || die ("Cannot open basic_tao.html file\n");

%mon2num = qw(
    January 01  February 02  March 03  April 04  May 05  June 06
    July 07  August 08  September 09  October 10 November 11 December 12
);

$date =~ /(.+) (.+), (.+) */;

$rev = "$3-$mon2num{$1}-$2";
print "Manual: tao-manual-$rev\n";

while (<F_IN>) {
  s/RRR/$rev/g;
  s/DDD/$date/g;
  print (F_OUT);
}

close (F_IN);
close (F_OUT);

`mv basic_tao.html /nfs/classe/www/html/bmad/`;
`cp tao.pdf /nfs/classe/www/html/bmad/tao-manual-$rev.pdf`;

