#!/usr/bin/perl

use File::Find;

#---------------------------------------------------------
# List of subroutines too low level to be mentioned

$tex_hash{"modulo2_sp"} = "modulo2_mod.f90";
$tex_hash{"modulo2_int"} = "modulo2_mod.f90";
$tex_hash{"modulo2_dp"} = "modulo2_mod.f90";
$tex_hash{"taylor_coef1"} = "bmad_taylor_mod.f90";
$tex_hash{"taylor_coef2"} = "bmad_taylor_mod.f90";
$tex_hash{"coord_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"amode_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"param_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"ran_uniform_vector"} = "random_mod.f90";
$tex_hash{"ran_uniform_scaler"} = "random_mod.f90";
$tex_hash{"integration_timer"} = "integration_timer_mod.f90";
$tex_hash{"test_f_em_field"} = "test_f_side.f90";
$tex_hash{"test_f_modes"} = "test_f_side.f90";
$tex_hash{"test_f_amode"} = "test_f_side.f90";
$tex_hash{"bmad_com_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"mode_info_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"lr_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr2_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr2_trans_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr1_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr1_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr2_long_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr2_long_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"lr_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"test_f_floor_position"} = "test_f_side.f90";
$tex_hash{"coord_in_orbit_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"modes_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"test_f_ele"} = "test_f_side.f90";
$tex_hash{"test_f_wig_term"} = "test_f_side.f90";
$tex_hash{"test_f_control"} = "test_f_side.f90";
$tex_hash{"test_f_sr1_wake"} = "test_f_side.f90";
$tex_hash{"test_f_coord"} = "test_f_side.f90";
$tex_hash{"test_f_taylor_term"} = "test_f_side.f90";
$tex_hash{"init_all_structs"} = "test_f_side.f90";
$tex_hash{"test_f_param"} = "test_f_side.f90";
$tex_hash{"test_f_ring"} = "test_f_side.f90";
$tex_hash{"test_f_mode_info"} = "test_f_side.f90";
$tex_hash{"test_f_wake"} = "test_f_side.f90";
$tex_hash{"test_f_lr_wake"} = "test_f_side.f90";
$tex_hash{"test_f_taylor"} = "test_f_side.f90";
$tex_hash{"test_f_bmad_com"} = "test_f_side.f90";
$tex_hash{"test_f_twiss"} = "test_f_side.f90";
$tex_hash{"test_f_linac_mode"} = "test_f_side.f90";
$tex_hash{"control_from_ring_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"taylor_term_in_taylor_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"wig_term_in_ele_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"ele_from_ring_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"taylor_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"floor_position_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"ele_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"linac_mode_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"ring_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"wig_term_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"twiss_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"em_field_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"taylor_term_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"orbit_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"control_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"parser_add_variable"} = "bmad_parser_mod.f90";
$tex_hash{"quad_mat2_calc"} = "make_mat6_mod.f90";
$tex_hash{"write_out"} = "io_mod.f90";
$tex_hash{"str"} = "io_mod.f90";
$tex_hash{"rchomp"} = "io_mod.f90";
$tex_hash{"qrfac"} = "lmdif_mod.f90";
$tex_hash{"if_error"} = "if_error.f90";
$tex_hash{"track_bend_edge"} = "track1_mod.f90";
$tex_hash{"ion_kick_2d"} = "ion_kick_2d.f90";
$tex_hash{"parser_add_lord"} = "bmad_parser_mod.f90";
$tex_hash{"mat6_multipole"} = "make_mat6_mod.f90";
$tex_hash{"mat4_multipole"} = "make_mat6_mod.f90";
$tex_hash{"evaluate_logical"} = "bmad_parser_mod.f90";
$tex_hash{"read_sr_wake"} = "bmad_parser_mod.f90";
$tex_hash{"read_lr_wake"} = "bmad_parser_mod.f90";
$tex_hash{"get_sequence_args"} = "bmad_parser_mod.f90";
$tex_hash{"allocate_pring"} = "bmad_parser_mod.f90";
$tex_hash{"add_taylor_term"} = "bmad_parser_mod.f90";
$tex_hash{"file_stack"} = "bmad_parser_mod.f90";
$tex_hash{"find_indexx"} = "bmad_parser_mod.f90";
$tex_hash{"get_attribute"} = "bmad_parser_mod.f90";
$tex_hash{"get_next_word"} = "bmad_parser_mod.f90";
$tex_hash{"get_overlay_group_names"} = "bmad_parser_mod.f90";
$tex_hash{"load_parse_line"} = "bmad_parser_mod.f90";
$tex_hash{"seq_expand1"} = "bmad_parser_mod.f90";
$tex_hash{"qromb_rad_int"} = "rad_int_common.f90";
$tex_hash{"derivs_bmad"} = "runge_kutta_mod.f90";
$tex_hash{"init_bmad_parser_common"} = "bmad_parser_mod.f90";
$tex_hash{"bbi_slice_calc"} = "make_mat6_mod.f90";
$tex_hash{"rkqs_bmad"} = "runge_kutta_mod.f90";
$tex_hash{"allocate_saved_orbit"} = "em_field_mod.f90";
$tex_hash{"preparse_element_init"} = "bmad_parser_mod.f90";
$tex_hash{"reallocate_control_"} = "bmad_utils_mod.f90";
$tex_hash{"compute_super_lord_s"} = "bmad_parser_mod.f90";
$tex_hash{"save_a_step"} = "em_field_mod.f90";
$tex_hash{"verify_valid_name"} = "bmad_parser_mod.f90";
$tex_hash{"warning"} = "bmad_parser_mod.f90";
$tex_hash{"bbi_kick_matrix"} = "make_mat6_bmad.f90";
$tex_hash{"track_solenoid_edge"} = "boris_mod.f90";
$tex_hash{"get_called_file"} = "bmad_parser_mod.f90";
$tex_hash{"mexp"} = "bmad_basic_mod.f90";
$tex_hash{"integration_timer_fibre"} = "integration_timer_mod.f90";
$tex_hash{"propagate_part_way"} = "rad_int_common.f90";
$tex_hash{"add_all_superimpose"} = "bmad_parser_mod.f90";
$tex_hash{"calc_g_params"} = "rad_int_common.f90";
$tex_hash{"diff"} = "integration_timer_mod.f90";
$tex_hash{"multi_turn_func"} = "multi_turn_tracking_to_mat.f90";
$tex_hash{"evaluate_value"} = "bmad_parser_mod.f90";
$tex_hash{"reallocate_bp_com_var"} = "bmad_parser_mod.f90";
$tex_hash{"error_exit"} = "bmad_parser_mod.f90";
$tex_hash{"do_vsp_eles"} = "create_vsp_volt_elements.f90";
$tex_hash{"word_to_value"} = "bmad_parser_mod.f90";
$tex_hash{"rkck_bmad"} = "runge_kutta_mod.f90";
$tex_hash{"term_diff"} = "integration_timer_mod.f90";
$tex_hash{"track1_bunch"} = "macro_particle_mod.f90";
$tex_hash{"pushit"} = "bmad_parser_mod.f90";
$tex_hash{"type_get"} = "bmad_parser_mod.f90";
$tex_hash{"get_taylor"} = "integration_timer_mod.f90";
$tex_hash{"order_macroparticles_in_z"} = "macro_particle_mod.f90";
$tex_hash{"integration_timer_ele"} = "integration_timer_mod.f90";
$tex_hash{"solenoid_mat_calc"} = "make_mat6_mod.f90";
$tex_hash{"lr_wake_add_to"} = "wake_mod.f90";
$tex_hash{"word_read"} =               "word_read.f90";
$tex_hash{"string_typ"} =              "string_typ.f90";
$tex_hash{"goodbye_exit"} =            "goodbye_exit.f90";
$tex_hash{"calc_file_number"} =        "calc_file_number.f90";
$tex_hash{"odd"} =                     "odd.f90";
$tex_hash{"fdjac2"} =                  "lmdif_mod.f90";
$tex_hash{"increment_file_number"} =   "increment_file_number.f90";
$tex_hash{"ion_kick"} =                "ion_kick.f90";
$tex_hash{"get_file_number"} =         "get_file_number.f90";
$tex_hash{"doubleup_quotes"} =         "doubleup_quotes.f90";
$tex_hash{"qrsolv"} =                  "lmdif_mod.f90";
$tex_hash{"string_to_int"} =           "string_to_int.f90";
$tex_hash{"match_word"} =              "match_word.f90";
$tex_hash{"ask"} =                     "ask.f90";
$tex_hash{"file_suffixer"} =           "file_suffixer.f90";
$tex_hash{"rel_typ"} =                 "rel_typ.f90";
$tex_hash{"change_file_number"} =      "change_file_number.f90";
$tex_hash{"settled_test"} =            "settled_test.f90";
$tex_hash{"spawn_command"} =           "spawn_command.f90";
$tex_hash{"to_node_name"} =            "to_node_name.f90";
$tex_hash{"file_get_open"} =           "file_get_open.f90";
$tex_hash{"match_word_exact"} =        "match_word_exact.f90";
$tex_hash{"file_directorizer"} =       "file_directorizer.f90";
$tex_hash{"make_legal_comment"} =      "make_legal_comment.f90";
$tex_hash{"ran_gauss_vector"} = "random_mod.f90";
$tex_hash{"ran_gauss_scaler"} = "random_mod.f90";
$tex_hash{"plot_it"} = "plot_example.f90";
$tex_hash{"compute_slave_aperture"} = "bookkeeper_mod.f90";


#---------------------------------------------------------
# make a list of names from bmad_subroutines.html

$tex_file =  '/home/dcs/new_lib/bmad/guide/subroutines.tex';
open (F_IN, $tex_file) || die ("Cannot open File: $tex_file");

while (<F_IN>) {
  if (/\\item\[(.*?)[\(\] ]/i)  {      # match to "\item[...]" 
    $name = $1;
    $name =~ s/protect\\parbox\{6in\}\{//;
    $name =~ tr/A-Z/a-z/;    # lowercase
    $name =~ s/\\//g;       # remove "\"
    chomp $name;
    if (! exists $tex_hash{$name}) { $tex_hash{$name} = $name; }
  }
}

## foreach $k (sort keys %tex_hash) {
##   print "$k  $tex_hash{$k}\n";
## }

#---------------------------------------------------------
# find all BMAD subroutines and see if they are in the tex file

print "\n---------------------------------------------------\n";
print "Subroutines in bmad but not in subroutines.tex:\n\n";

find (\&searchit, '/home/dcs/new_lib/bmad');
find (\&searchit, '/home/dcs/new_lib/dcslib');
find (\&searchit, '/home/dcs/new_lib/cesr_utils');

#---------------------------------------------------------
# find all BMAD subroutines and see if they are in the tex file

print "\n---------------------------------------------------\n";
print "Subroutines in bmad but not in subroutines.tex:\n\n";

while (($key, $val) = each %f90_hash) {
  if (! exists $tex_hash{$key}) {
#    printf "%-30s%s\n", ":$key:", $val;
    printf "\$tex_hash\{\"$key\"\} = \"$val\"\;\n";
  }
}

#---------------------------------------------------------
# find all listings in the tex file and see if there are bmad subroutines

print "\n---------------------------------------------------\n";
print "Subroutines in subroutines.tex or bmad_subroutines.pl that are not in bmad:\n\n";

while (($key, $val) = each %tex_hash) {
  if (! exists $f90_hash{$key}) {
    print ":$key:\n";
  }
}

#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------

sub searchit {

  if (!/\.f90$/) {return;}

  $file = $_;

  open (F_IN, $file) || die ("Cannot open File: $_");

  while (<F_IN>) {
    $now = $_;

    if (/^ *interface *$/i) {    # skip interface blocks
      while (<F_IN>) {
        if (/^ *end interface/i) {last;}
      }
    }

# Look for "! subroutine ..." after "!+" to get description.

    if ($old =~ /^\!\+/) {
      if (/^\! *subroutine */i || /^\! *function */i) {
      
        $this = $';     # strip off "subroutine"
        chomp $this;
        $this =~ tr/A-Z/a-z/; #lowercase
        $name = $this;
        $name =~ s/ *\(.*//;    # strip off " (..."

        if (exists $tex_hash{$name}) {
##          print "$this  In list \n";
        } else {
          print "\nFile: $file\n";
          $this2 = $this; 
          $this2 =~ s/\s*\(.*//;
          print "\\index\{Routine\!$this2\}\n";
          $this =~ s/_/\\_/g;      # "_"   -->  "\_"
          print "\\item\[$this\] \\Newline \n";
          $l2 = <F_IN>; 
          $l2 =~ s/_/\\_/g;      # "_"   -->  "\_"
          if (! ($l2 =~ /^\! *[\n]$/)) {$l2 =~ s/^\! //; print "$l2"; $temp = <F_IN>;}
          $l2 = <F_IN>; 
          $l2 =~ s/_/\\_/g;      # "_"   -->  "\_"
          if (! ($l2 =~ /^\! *[\n]$/)) {$l2 =~ s/^\! //; print "$l2";}
          $l2 = <F_IN>; 
          $l2 =~ s/_/\\_/g;      # "_"   -->  "\_"
          if (! ($l2 =~ /^\! *[\n]$/)) {$l2 =~ s/^\! //; print "$l2";}
        }
      }
    }

# Look for "subroutine ..." to get name list.

    if (/^ *subroutine /i || /^ *recursive subroutine /i || 
        /^ *function /i   || /^ *real\(rp\) *function /i || 
        /^ *elemental subroutine /i || /^ *interface /i) {
      $name = $';              # strip off "subroutine
      $name =~ s/ *\(.*//;     # strip off " (..."
      $name =~ s/ +$//;        # strip off trailing blank if no arg list
      $name =~ tr/A-Z/a-z/; #lowercase
      chomp $name;
      if (! exists $f90_hash{$name}) { $f90_hash{$name} = $file; }

      # skip rest of routine including contained routines

      $count = 1;
      while (<F_IN>) {
        if (/^ *end /i) {
          $_ = $';  
          if (/^ *subroutine/i || /^ *function/i || /^ *interface/i) {
            $count = $count - 1;
          }
        }
        elsif (/^ *subroutine /i || /^ *recursive subroutine /i || 
              /^ *function /i || /^ *elemental subroutine /i ||
              /^ *real\(rp\) *function /i || /^ *interface /i) {
          $count = $count + 1;
        }
        if ($count == 0) {last;}
      }

    }

    $old = $_;
  }

  close (F_IN);

  $_ = $file;

}
