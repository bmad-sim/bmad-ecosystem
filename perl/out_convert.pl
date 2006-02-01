#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/map_no_offsets/map_with_offsets/g;
    s/MAP_NO_OFFSETS/MAP_WITH_OFFSETS/g;
    s/exact_rad_int_calc/map_with_offsets/g;

    {print (F_OUT);}

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
