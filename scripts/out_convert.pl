#!/usr/bin/perl   

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/control_type \=\= overlay_lord\$/lord_status \=\= overlay_lord\$/g;
    s/control_type \/\= overlay_lord\$/lord_status \/\= overlay_lord\$/g;
    s/control_type \= overlay_lord\$/lord_status \= overlay_lord\$/g;

    s/control_type \=\= super_lord\$/lord_status \=\= super_lord\$/g;
    s/control_type \/\= super_lord\$/lord_status \/\= super_lord\$/g;
    s/control_type \= super_lord\$/lord_status \= super_lord\$/g;

    s/control_type \=\= multipass_lord\$/lord_status \=\= multipass_lord\$/g;
    s/control_type \/\= multipass_lord\$/lord_status \/\= multipass_lord\$/g;
    s/control_type \= multipass_lord\$/lord_status \= multipass_lord\$/g;

    s/control_type \=\= girder_lord\$/lord_status \=\= girder_lord\$/g;
    s/control_type \/\= girder_lord\$/lord_status \/\= girder_lord\$/g;
    s/control_type \= girder_lord\$/lord_status \= girder_lord\$/g;

    s/control_type \=\= group_lord\$/lord_status \=\= group_lord\$/g;
    s/control_type \/\= group_lord\$/lord_status \/\= group_lord\$/g;
    s/control_type \= group_lord\$/lord_status \= group_lord\$/g;

    s/control_type \=\= overlay_slave\$/slave_status \=\= overlay_slave\$/g;
    s/control_type \/\= overlay_slave\$/slave_status \/\= overlay_slave\$/g;
    s/control_type \= overlay_slave\$/slave_status \= overlay_slave\$/g;

    s/control_type \=\= super_slave\$/slave_status \=\= super_slave\$/g;
    s/control_type \/\= super_slave\$/slave_status \/\= super_slave\$/g;
    s/control_type \= super_slave\$/slave_status \= super_slave\$/g;

    s/control_type \=\= multipass_slave\$/slave_status \=\= multipass_slave\$/g;
    s/control_type \/\= multipass_slave\$/slave_status \/\= multipass_slave\$/g;
    s/control_type \= multipass_slave\$/slave_status \= multipass_slave\$/g;

    print (F_OUT);
  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
