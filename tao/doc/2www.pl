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

`python3 generate_tex_for_pipe_commands.py`;
`pdflatex tao`;

$file = "tao_template.html";
open (F_IN, $file) || die ("Cannot open File: $file\n");

open (F_OUT, ">basic_tao.html") || die ("Cannot open basic_tao.html file\n");

%month2num = qw(
    January 01  February 02  March 03  April 04  May 05  June 06
    July 07  August 08  September 09  October 10 November 11 December 12
);

$date =~ /(.+) (.+), (.+) */;
$day = $2;
if ($day < 10) { $day = "0$day"; }    # Two digit day

$rev = "$3-$month2num{$1}-$day";
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
`python3 generate_tex_for_pipe_commands.py`
