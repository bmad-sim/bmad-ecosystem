#!/usr/bin/perl   

use File::Copy;

foreach $file (@ARGV) {  
  print "Converting: $file\n";
  copy( $file, "$file~" ) or warn "Couldn't copy $file to $file~:\n $!\n"; 
  print "Copied current version to $file~\n";
  open (F_IN, $file) || die ("Cannot Open File");
  open (F_OUT, ">temp.out");

  while (<F_IN>) {

    s/\(macro_init_twiss_struct\) x, y/\(macro_init_twiss_struct\) a, b/;

    s/ring_struct/lat_struct/ig;
    s/\%control_([\;\(\,\)\s])/\%control\1/ig;
    s/\%ic_([\;\(\,\)\s])/\%ic\1/ig;
    s/\%ele_([\;\(\,\)\s\%])/\%ele\1/ig;
    s/\%n_ele_ring/\%n_ele_use/ig;
    s/\%n_ele_use/\%n_ele_track/ig;
    s/\%n_ele_lat/\%n_ele_track/ig;
    s/floor\%a/floor\%x/ig;
    s/floor\%b/floor\%y/ig;
    s/BEAM_ENERGY/E_TOT/ig;
    s/beam_energy/E_tot/ig;
    s/ENERGY_START/E_TOT_START/ig;
    s/energy_start/E_tot_start/ig;
    s/sr2_long/sr_mode_long/ig;
    s/sr2_trans/sr_mode_trans/ig;
    s/sr1/sr_table/ig;
    s/sr2/sr_mode/ig;
    s/ix_ring/ix_lat/ig;
    s/field\%b_pole/field\%B/ig;
    ## s/key/class/ig;

    s/sr1_wake_struct/sr_table_wake_struct/ig;
    s/sr2_wake_struct/sr_mode_wake_struct/ig;

    s/bmad_com_struct/bmad_common_struct/ig;
    s/bp_com_struct/bp_common_struct/ig;
    s/ptc_com_struct/ptc_common_struct/ig;
    s/synch_rad_com_struct/synch_rad_common_struct/ig;
    s/modes_struct/normal_modes_struct/ig;
    s/mode_struct/normal_mode_struct/ig;
    s/param_struct/lat_param_struct/ig;

    s/transfer_ring_taylors/transfer_lat_taylors/ig;
    s/test_f_ring/test_f_lat/ig;
    s/parser_ring_struct/parser_lat_struct/ig;
    s/transfer_ring_parameters/transfer_lat_parameters/ig;
    s/init_ring/init_lat/ig;
    s/clear_ring_1turn_mats/clear_lat_1turn_mats/ig;
    s/transfer_ring/transfer_lat/ig;
    s/allocate_ring_ele_/allocate_lat_ele_/ig;
    s/deallocate_ring_pointers/deallocate_lat_pointers/ig;
    s/transfer_ring_taylors/transfer_lat_taylors/ig;
    s/ring_equal_ring/lat_equal_ring/ig;
    s/ring_vec_equal_ring_vec/lat_vec_equal_ring_vec/ig;
    s/ring_to_layout/lat_to_layout/ig;
    s/ring_reverse/lat_reverse/ig;
    s/C_ring/C_lat/ig;
    s/ring_to_c/lat_to_c/ig;
    s/ring_to_f2/lat_to_f2/ig;
    s/ele_from_ring_to_f2/ele_from_lat_to_f2/ig;
    s/control_from_ring_to_f2/control_from_lat_to_f2/ig;
    s/ring_geometry/lat_geometry/ig;
    s/ring_make_mat6/lat_make_mat6/ig;
    s/check_ring_controls/check_lat_controls/ig;
    s/split_ring/split_lat/ig;
    s/compress_ring/compress_lat/ig;
    s/make_hybrid_ring/make_hybrid_lat/ig;
    s/ring_to_quad_calib/lat_to_quad_calib/ig;
    s/ring_master_struct/lat_master_struct/ig;

    s/beta_x/beta_a/g;
    s/Beta_x/Beta_a/g;
    s/BETA_X/BETA_A/g;

    s/beta_y/beta_b/g;
    s/Beta_y/Beta_b/g;
    s/BETA_Y/BETA_B/g;

    s/alpha_x/alpha_a/g;
    s/Alpha_x/Alpha_a/g;
    s/ALPHA_X/ALPHA_A/g;

    s/alpha_y/alpha_b/g;
    s/Alpha_y/Alpha_b/g;
    s/ALPHA_Y/ALPHA_B/g;

    if (/theta/i) {print (F_OUT); next;}

    s/eta_x/eta_a/g;
    s/Eta_x/Eta_a/g;
    s/ETA_X/ETA_A/g;

    s/eta_y/eta_b/g;
    s/Eta_y/Eta_b/g;
    s/ETA_Y/ETA_B/g;

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

    if (/mag\%wi\%w/) {print (F_OUT); next;}
    if (/aperture\%/) {print (F_OUT); next;}
    if (/\%x\%min/) {print (F_OUT); next;}
    if (/\%y\%min/) {print (F_OUT); next;}
    if (/\%x\%max/) {print (F_OUT); next;}
    if (/\%y\%max/) {print (F_OUT); next;}
    if (/\%x\%draw/) {print (F_OUT); next;}
    if (/\%y\%draw/) {print (F_OUT); next;}
    if (/\%x\%major/) {print (F_OUT); next;}
    if (/\%y\%major/) {print (F_OUT); next;}
    if (/\%x\%places/) {print (F_OUT); next;}
    if (/\%y\%places/) {print (F_OUT); next;}
    if (/origin\%x/) {print (F_OUT); next;}
    if (/origin\%y/) {print (F_OUT); next;}
    if (/plot\%x/) {print (F_OUT); next;}
    if (/plot\%y/) {print (F_OUT); next;}
    if (/floor\%x/) {print (F_OUT); next;}
    if (/floor\%y/) {print (F_OUT); next;}
    if (/graph\%x/) {print (F_OUT); next;}
    if (/graph\%y/) {print (F_OUT); next;}
    if (/plot\(.*\)\%x/) {print (F_OUT); next;}
    if (/plot\(.*\)\%y/) {print (F_OUT); next;}
    if (/graph\(.*\)\%x/) {print (F_OUT); next;}
    if (/graph\(.*\)\%y/) {print (F_OUT); next;}
    if (/title\(.*\)\%x/) {print (F_OUT); next;}
    if (/title\(.*\)\%y/) {print (F_OUT); next;}

    s/\%a([\;\(\,\)\s])/\%a_pole\1/ig;
    s/\%b([\;\(\,\)\s])/\%b_pole\1/ig;
    s/\%x([\;\%\s\)\,])/\%a\1/ig;
    s/\%y([\;\%\s\)\,])/\%b\1/ig;
    s/a\%eta_lab/x\%eta/ig;
    s/b\%eta_lab/y\%eta/ig;
    s/a\%etap_lab/x\%etap/ig;
    s/b\%etap_lab/y\%etap/ig;    

    print (F_OUT);

  }

  close (F_IN);
  close (F_OUT);

  $file2 = $file;
  move( "temp.out", $file2 ) or warn "Couldn't move temp.out to $file2: $!\n";

}
