#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/plot_page\%region/plot_region/g;

    print (F_OUT);

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
