#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/E-05/E-08/g;
    s/E-04/E-07/g;
    s/E-03/E-06/g;
    s/E-02/E-05/g;
    s/E-01/E-04/g;
    s/E\+00/E-03/g;

    print (F_OUT);

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
