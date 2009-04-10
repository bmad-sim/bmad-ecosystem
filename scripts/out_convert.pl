#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/sr_interface/synrad_interface/;
    s/sr_struct/synrad_struct/;
    s/sr_mod/synrad_mod/;
    s/calculate_sr_power/calculate_synrad_power/;

    print (F_OUT);
  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
