#!/usr/bin/perl   

$n = @ARGV;
if ($n == 0) {
  print "Usage:\n";
  print "  input_file_convert.pl  <tune_scan input file>\n";
  exit;
}

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    if (!/Q_z0/i && !/Q_z1/i && !/dQ_z/) {
      s/Q_z/scan_params\%Q_z/i;
    }
    s/lat_file/scan_params\%lat_file/i;
    s/n_turn/scan_params\%n_turn/i;
    s/n_part/scan_params\%n_part/i;
    s/init/scan_params\%init/i;
    s/particle/scan_params\%particle/i;
    s/i_train/scan_params\%i_train/i;
    s/j_car/scan_params\%j_car/i;
    s/n_trains_tot/scan_params\%n_trains_tot/i;
    s/n_cars/scan_params\%n_cars/i;
    s/current/scan_params\%current/i;
    s/lrbbi/scan_params\%lrbbi/i;
    s/beambeam_ip/scan_params\%beambeam_ip/i;
    s/close_pretz/scan_params\%close_pretz/i;
    s/close_vert/scan_params\%close_vert/i;
    s/rec_taylor/scan_params\%rec_taylor/i;
    s/slices/scan_params\%slices/i;

    print (F_OUT);

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  rename "temp.out", $file2;

}
