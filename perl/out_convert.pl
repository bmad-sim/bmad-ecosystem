#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/etap_x/etap_a/g;
    s/Etap_x/Etap_a/g;
    s/ETAP_X/ETAP_A/g;

    s/etap_y/etap_b/g;
    s/Etap_y/Etap_b/g;
    s/ETAP_Y/ETAP_B/g;

    s/phi_x/phi_a/g;
    s/Phi_x/Phi_a/g;
    s/PHI_X/PHI_A/g;

    s/phi_y/phi_b/g;
    s/Phi_y/Phi_b/g;
    s/PHI_Y/PHI_B/g;

    {print (F_OUT);}

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
