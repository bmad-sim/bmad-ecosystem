#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    ## if (/calculation/) {print "Calc found: $_";}
    s/a\%eta_lab/x\%eta/g;
    s/a\%etap_lab/x\%etap/g;
    s/b\%eta_lab/y\%eta/g;
    s/b\%etap_lab/y\%etap/g;

    print (F_OUT);

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
