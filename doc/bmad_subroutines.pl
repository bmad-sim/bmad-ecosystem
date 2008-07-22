#!/usr/bin/perl

use File::Find;

#---------------------------------------------------------
# List of subroutines too low level to be mentioned

$pl_hash{"lsc_y0_kick_calc"} = "csr_mod.f90";
$pl_hash{"test_f_xy_disp"} = "test_f_side.f90";
$pl_hash{"re_allocate2_real"} = "re_allocate_mod.f90";
$pl_hash{"this_bookkeeper"} = "bookkeeper_mod.f90";
$pl_hash{"re_allocate2_integer"} = "re_allocate_mod.f90";
$pl_hash{"xy_disp_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"re_allocate2_string"} = "re_allocate_mod.f90";
$pl_hash{"set_string"} = "fortran_and_cpp.f90";
$pl_hash{"re_allocate2_logical"} = "re_allocate_mod.f90";
$pl_hash{"form_digested_bmad_file_name"} = "bmad_parser_mod.f90";
$pl_hash{"reuse_taylor_elements"} = "bmad_parser_mod.f90";
$pl_hash{"add_this_multipass"} = "bmad_parser_mod.f90";
$pl_hash{"settable_dep_var_bookkeeping"} = "bmad_parser_mod.f90";
$pl_hash{"find_this_file"} = "bmad_parser_mod.f90";
$pl_hash{"parser_set_ele_defaults"} = "bmad_parser_mod.f90";
$pl_hash{"save_taylor_elements"} = "bmad_parser_mod.f90";
$pl_hash{"initialize_pauli_vector"} = "spin_mod.f90";
$pl_hash{"track1_macro_sr_trans_wake"} = "macroparticle_mod.f90";
$pl_hash{"track1_macro_lr_wake"} = "macroparticle_mod.f90";
$pl_hash{"makeup_overlay_and_girder_slave"} = "bookkeeper_mod.f90";
$pl_hash{"compute_slave_coupler"} = "bookkeeper_mod.f90";
$pl_hash{"makeup_multipass_slave"} = "bookkeeper_mod.f90";
$pl_hash{"modulo2_sp"} = "modulo2_mod.f90";
$pl_hash{"modulo2_int"} = "modulo2_mod.f90";
$pl_hash{"modulo2_dp"} = "modulo2_mod.f90";
$pl_hash{"taylor_coef1"} = "bmad_taylor_mod.f90";
$pl_hash{"taylor_coef2"} = "bmad_taylor_mod.f90";
$pl_hash{"coord_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"amode_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"param_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"ran_uniform_vector"} = "random_mod.f90";
$pl_hash{"ran_uniform_scaler"} = "random_mod.f90";
$pl_hash{"integration_timer"} = "integration_timer_mod.f90";
$pl_hash{"test_f_em_field"} = "test_f_side.f90";
$pl_hash{"test_f_modes"} = "test_f_side.f90";
$pl_hash{"test_f_amode"} = "test_f_side.f90";
$pl_hash{"bmad_com_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"mode_info_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"lr_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"sr_mode_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"sr_mode_trans_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"sr_table_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"sr_table_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"sr_mode_long_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"sr_mode_long_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"lr_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"test_f_floor_position"} = "test_f_side.f90";
$pl_hash{"modes_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"test_f_ele"} = "test_f_side.f90";
$pl_hash{"test_f_wig_term"} = "test_f_side.f90";
$pl_hash{"test_f_control"} = "test_f_side.f90";
$pl_hash{"test_f_sr_table_wake"} = "test_f_side.f90";
$pl_hash{"test_f_coord"} = "test_f_side.f90";
$pl_hash{"test_f_taylor_term"} = "test_f_side.f90";
$pl_hash{"init_all_structs"} = "test_f_side.f90";
$pl_hash{"test_f_param"} = "test_f_side.f90";
$pl_hash{"test_f_lat"} = "test_f_side.f90";
$pl_hash{"test_f_mode_info"} = "test_f_side.f90";
$pl_hash{"test_f_wake"} = "test_f_side.f90";
$pl_hash{"test_f_lr_wake"} = "test_f_side.f90";
$pl_hash{"test_f_taylor"} = "test_f_side.f90";
$pl_hash{"test_f_bmad_com"} = "test_f_side.f90";
$pl_hash{"test_f_twiss"} = "test_f_side.f90";
$pl_hash{"test_f_linac_mode"} = "test_f_side.f90";
$pl_hash{"control_from_lat_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"taylor_term_in_taylor_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"wig_term_in_ele_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"ele_from_lat_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"taylor_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"floor_position_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"ele_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"linac_mode_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"lat_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"wake_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"wig_term_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"twiss_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"em_field_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"taylor_term_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"control_to_f2"} = "bmad_and_cpp.f90";
$pl_hash{"parser_add_variable"} = "bmad_parser_mod.f90";
$pl_hash{"quad_mat2_calc"} = "make_mat6_mod.f90";
$pl_hash{"write_out"} = "io_mod.f90";
$pl_hash{"str"} = "io_mod.f90";
$pl_hash{"rchomp"} = "io_mod.f90";
$pl_hash{"qrfac"} = "lmdif_mod.f90";
$pl_hash{"if_error"} = "if_error.f90";
$pl_hash{"track_bend_edge"} = "track1_mod.f90";
$pl_hash{"ion_kick_2d"} = "ion_kick_2d.f90";
$pl_hash{"parser_add_lord"} = "bmad_parser_mod.f90";
$pl_hash{"mat6_multipole"} = "make_mat6_mod.f90";
$pl_hash{"mat4_multipole"} = "make_mat6_mod.f90";
$pl_hash{"evaluate_logical"} = "bmad_parser_mod.f90";
$pl_hash{"read_sr_wake"} = "bmad_parser_mod.f90";
$pl_hash{"read_lr_wake"} = "bmad_parser_mod.f90";
$pl_hash{"get_sequence_args"} = "bmad_parser_mod.f90";
$pl_hash{"allocate_plat"} = "bmad_parser_mod.f90";
$pl_hash{"add_taylor_term"} = "bmad_parser_mod.f90";
$pl_hash{"file_stack"} = "bmad_parser_mod.f90";
$pl_hash{"find_indexx"} = "bmad_parser_mod.f90";
$pl_hash{"get_attribute"} = "bmad_parser_mod.f90";
$pl_hash{"get_next_word"} = "bmad_parser_mod.f90";
$pl_hash{"get_overlay_group_names"} = "bmad_parser_mod.f90";
$pl_hash{"load_parse_line"} = "bmad_parser_mod.f90";
$pl_hash{"seq_expand1"} = "bmad_parser_mod.f90";
$pl_hash{"qromb_rad_int"} = "rad_int_common.f90";
$pl_hash{"init_bmad_parser_common"} = "bmad_parser_mod.f90";
$pl_hash{"bbi_slice_calc"} = "make_mat6_mod.f90";
$pl_hash{"rkqs_bmad"} = "runge_kutta_mod.f90";
$pl_hash{"allocate_saved_orbit"} = "em_field_mod.f90";
$pl_hash{"preparse_element_init"} = "bmad_parser_mod.f90";
$pl_hash{"reallocate_control"} = "bmad_utils_mod.f90";
$pl_hash{"compute_super_lord_s"} = "bmad_parser_mod.f90";
$pl_hash{"save_a_step"} = "em_field_mod.f90";
$pl_hash{"verify_valid_name"} = "bmad_parser_mod.f90";
$pl_hash{"warning"} = "bmad_parser_mod.f90";
$pl_hash{"bbi_kick_matrix"} = "make_mat6_bmad.f90";
$pl_hash{"track_solenoid_edge"} = "boris_mod.f90";
$pl_hash{"get_called_file"} = "bmad_parser_mod.f90";
$pl_hash{"mexp"} = "bmad_basic_mod.f90";
$pl_hash{"integration_timer_fibre"} = "integration_timer_mod.f90";
$pl_hash{"propagate_part_way"} = "rad_int_common.f90";
$pl_hash{"add_all_superimpose"} = "bmad_parser_mod.f90";
$pl_hash{"diff"} = "integration_timer_mod.f90";
$pl_hash{"multi_turn_func"} = "multi_turn_tracking_to_mat.f90";
$pl_hash{"evaluate_value"} = "bmad_parser_mod.f90";
$pl_hash{"reallocate_bp_com_var"} = "bmad_parser_mod.f90";
$pl_hash{"parser_debug_print_info"} = "bmad_parser_mod.f90";
$pl_hash{"parser_expand_line"} = "bmad_parser_mod.f90";
$pl_hash{"error_exit"} = "bmad_parser_mod.f90";
$pl_hash{"do_vsp_eles"} = "create_vsp_volt_elements.f90";
$pl_hash{"word_to_value"} = "bmad_parser_mod.f90";
$pl_hash{"rkck_bmad"} = "runge_kutta_mod.f90";
$pl_hash{"term_diff"} = "integration_timer_mod.f90";
$pl_hash{"pushit"} = "bmad_parser_mod.f90";
$pl_hash{"type_get"} = "bmad_parser_mod.f90";
$pl_hash{"get_taylor"} = "integration_timer_mod.f90";
$pl_hash{"order_macroparticles_in_z"} = "macro_particle_mod.f90";
$pl_hash{"integration_timer_ele"} = "integration_timer_mod.f90";
$pl_hash{"solenoid_mat_calc"} = "make_mat6_mod.f90";
$pl_hash{"lr_wake_add_to"} = "wake_mod.f90";
$pl_hash{"word_read"} =               "word_read.f90";
$pl_hash{"string_typ"} =              "string_typ.f90";
$pl_hash{"goodbye_exit"} =            "goodbye_exit.f90";
$pl_hash{"calc_file_number"} =        "calc_file_number.f90";
$pl_hash{"odd"} =                     "odd.f90";
$pl_hash{"fdjac2"} =                  "lmdif_mod.f90";
$pl_hash{"increment_file_number"} =   "increment_file_number.f90";
$pl_hash{"ion_kick"} =                "ion_kick.f90";
$pl_hash{"get_file_number"} =         "get_file_number.f90";
$pl_hash{"doubleup_quotes"} =         "doubleup_quotes.f90";
$pl_hash{"qrsolv"} =                  "lmdif_mod.f90";
$pl_hash{"match_word"} =              "match_word.f90";
$pl_hash{"ask"} =                     "ask.f90";
$pl_hash{"file_suffixer"} =           "file_suffixer.f90";
$pl_hash{"rel_typ"} =                 "rel_typ.f90";
$pl_hash{"change_file_number"} =      "change_file_number.f90";
$pl_hash{"settled_test"} =            "settled_test.f90";
$pl_hash{"spawn_command"} =           "spawn_command.f90";
$pl_hash{"to_node_name"} =            "to_node_name.f90";
$pl_hash{"file_get_open"} =           "file_get_open.f90";
$pl_hash{"file_directorizer"} =       "file_directorizer.f90";
$pl_hash{"make_legal_comment"} =      "make_legal_comment.f90";
$pl_hash{"ran_gauss_vector"} = "random_mod.f90";
$pl_hash{"ran_gauss_scaler"} = "random_mod.f90";
$pl_hash{"plot_it"} = "plot_example.f90";
$pl_hash{"compute_slave_aperture"} = "bookkeeper_mod.f90";
$pl_hash{"delete_lattice_control_struct"} = "delete_lattice_control_struct.f90";
$pl_hash{"test_f_sr_mode_wake"} = "test_f_side.f90";
$pl_hash{"norm66"} = "eigen_mod.f90";
$pl_hash{"find_format"} = "output_mod.f90";
$pl_hash{"icm_typ"} = "icm_typ.f90";
$pl_hash{"bjmt_int_v"} = "ibs_mod.f90";
$pl_hash{"probability_funct"} = "probability_funct.f90";
$pl_hash{"slice_ele_calc"} = "slice_ele_calc.f90";
$pl_hash{"bmad_to_mad_or_xsif"} = "io_mod.f90";
$pl_hash{"integer_read"} = "integer_read.f90";
$pl_hash{"output_real"} = "output_mod.f90";
$pl_hash{"re_associate_real"} = "reallocate_mod.f90";
$pl_hash{"run_timer"} = "run_timer.f90";
$pl_hash{"str_upcase"} = "str_upcase.f90";
$pl_hash{"mat6_add_pitch"} = "bmad_utils_mod.f90";
$pl_hash{"boris_energy_correction"} = "boris_mod.f90";
$pl_hash{"re_associate_string"} = "reallocate_mod.f90";
$pl_hash{"sbend_body_with_k1_map"} = "make_mat6_mod.f90";
$pl_hash{"fullfilename"} = "fullfilename.f90";
$pl_hash{"bjmt_int_p"} = "ibs_mod.f90";
$pl_hash{"output_int"} = "output_mod.f90";
$pl_hash{"transfer_map_calc"} = "transfer_map_calc.f90";
$pl_hash{"re_allocate_real"} = "reallocate_mod.f90";
$pl_hash{"transfer_wake"} = "bmad_utils_mod.f90";
$pl_hash{"ordersys"} = "eigen_mod.f90";
$pl_hash{"i_size"} = "bmad_and_cpp.f90";
$pl_hash{"indexx"} = "indexx_mod.f90";
$pl_hash{"output_lines"} = "output_mod.f90";
$pl_hash{"lmpar"} = "lmdif_mod.f90";
$pl_hash{"create_nir_shuntcur_elements"} = "cesr_basic_mod.f90";
$pl_hash{"brent2"} = "brent2.f90";
$pl_hash{"ety"} = "eigen_mod.f90";
$pl_hash{"save_bunch_track"} = "save_bunch_track.f90";
$pl_hash{"normalized_quaternion"} = "spin_mod.f90";
$pl_hash{"bjmt_int_h"} = "ibs_mod.f90";
$pl_hash{"value_to_line"} = "io_mod.f90";
$pl_hash{"spin_omega_at"} = "spin_mod.f90";
$pl_hash{"eigensys"} = "eigen_mod.f90";
$pl_hash{"re_allocate_string"} = "reallocate_mod.f90";
$pl_hash{"query_string"} = "query_string.f90";
$pl_hash{"output_direct"} = "output_mod.f90";
$pl_hash{"header_io"} = "output_mod.f90";
$pl_hash{"output_line4"} = "output_mod.f90";
$pl_hash{"element_out"} = "io_mod.f90";
$pl_hash{"output_logical"} = "output_mod.f90";
$pl_hash{"r_size"} = "bmad_and_cpp.f90";
$pl_hash{"track1_trans_space_charge"} = "trans_space_charge_mod.f90";
$pl_hash{"operator"} = "equality_mod.f90";
$pl_hash{"calc_superimpose_key"} = "add_superimpose.f90";
$pl_hash{"integrand"} = "ibs_mod.f90";
$pl_hash{"make_mat6_trans_space_charge"} = "trans_space_charge_mod.f90";
$pl_hash{"exp_bessi0"} = "touschek_mod.f90";
$pl_hash{"re_associate_logical"} = "reallocate_mod.f90";
$pl_hash{"transfer_rad_int_struct"} = "rad_int_common.f90";
$pl_hash{"eval_logical"} = "utilities_mod.f90";
$pl_hash{"file_get"} = "file_get.f90";
$pl_hash{"real_read"} = "real_read.f90";
$pl_hash{"lmdif"} = "lmdif_mod.f90";
$pl_hash{"grad_loss_macro_sr_wake_calc"} = "macroparticle_mod.f90";
$pl_hash{"sr_table_add_long_kick"} = "wake_mod.f90";
$pl_hash{"to_word"} = "to_word.f90";
$pl_hash{"csr_bell"} = "csr_bell.f90";
$pl_hash{"word_len"} = "word_len.f90";
$pl_hash{"add_lattice_control_structs"} = "add_lattice_control_structs.f90";
$pl_hash{"ety2"} = "eigen_mod.f90";
$pl_hash{"init_coord_struct"} = "bmad_and_cpp.f90";
$pl_hash{"etyt"} = "eigen_mod.f90";
$pl_hash{"out_io_line_out"} = "output_mod.f90";
$pl_hash{"etdiv"} = "eigen_mod.f90";
$pl_hash{"mp_slice_equal_mp_slice"} = "macroparticle_mod.f90";
$pl_hash{"insert_numbers"} = "output_mod.f90";
$pl_hash{"int_typ"} = "int_typ.f90";
$pl_hash{"print_eq_ele"} = "equality_mod.f90";
$pl_hash{"inverse"} = "inverse.f90";
$pl_hash{"track1_macro_bunch"} = "macroparticle_mod.f90";
$pl_hash{"re_allocate_integer"} = "reallocate_mod.f90";
$pl_hash{"calc_wiggler_g_params"} = "rad_int_common.f90";
$pl_hash{"find1_indexx"} = "io_mod.f90";
$pl_hash{"mp_bunch_equal_mp_bunch"} = "macroparticle_mod.f90";
$pl_hash{"ele_geometry"} = "lat_geometry.f90";
$pl_hash{"re_allocate_logical"} = "reallocate_mod.f90";
$pl_hash{"mat6_add_offsets"} = "mat6_add_offsets.f90";
$pl_hash{"find_file"} = "find_file.f90";
$pl_hash{"assignment"} = "bmad_taylor_mod.f90";
$pl_hash{"re_associate_integer"} = "reallocate_mod.f90";
$pl_hash{"logical_typ"} = "logical_typ.f90";

#---------------------------------------------------------
# make a list of names from bmad_subroutines.html

$tex_file =  'subroutines.tex';
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

find (\&searchit, '..');
find (\&searchit, '../../dcslib');
find (\&searchit, '../../cesr_utils');

#---------------------------------------------------------
# find all BMAD subroutines and see if they are in the tex file

print "\n---------------------------------------------------\n";
print "Subroutines in bmad but not in subroutines.tex:\n\n";

while (($key, $val) = each %f90_hash) {
  if (! exists $tex_hash{$key} && ! exists $pl_hash{$key}) {
#    printf "%-30s%s\n", ":$key:", $val;
    printf "\$pl_hash\{\"$key\"\} = \"$val\"\;\n";
  }
}

#---------------------------------------------------------
# find all listings in the tex file and see if there are bmad subroutines

print "\n---------------------------------------------------\n";
print "Subroutines in subroutines.tex that are not in bmad:\n\n";

while (($key, $val) = each %tex_hash) {
  if (! exists $f90_hash{$key}) {
    print ":$key:\n";
  }
}

#---------------------------------------------------------
# find all listings in the tex file and see if there are bmad subroutines

print "\n---------------------------------------------------\n";
print "Subroutines in bmad_subroutines.pl that are not in bmad:\n\n";

while (($key, $val) = each %pl_hash) {
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
##  print "File: $file\n";

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
      
        $this = $';     #' strip off "subroutine"
        chomp $this;
        $this =~ tr/A-Z/a-z/; #lowercase
        $name = $this;
        $name =~ s/ *\(.*//;    # strip off " (..."
        $name =~ s/ +$//;       # strip off trailing blanks if no arg list.
        # print "Found: $name\n";
        if (!exists $tex_hash{$name} && !exists $pl_hash{$name}) {
          print "\nFile: $file\n";
          $this2 = $this; 
          $this2 =~ s/\s*\(.*//;
          print "\\index\{Routine\!$this2\}\n";
          print "\\item\[$this\] \\Newline \n";
          $l2 = <F_IN>; 
          $l2 =~ s/\%/\\\%/g;
          if (! ($l2 =~ /^\! *[\n]$/)) {$l2 =~ s/^\! //; print "$l2"; $temp = <F_IN>;}
          $l2 = <F_IN>; 
          $l2 =~ s/\%/\\\%/g;
          if (! ($l2 =~ /^\! *[\n]$/)) {$l2 =~ s/^\! //; print "$l2";}
          $l2 = <F_IN>; 
          $l2 =~ s/\%/\\\%/g;
          if (! ($l2 =~ /^\! *[\n]$/)) {$l2 =~ s/^\! //; print "$l2";}
        }
      }
    }

# Look for "subroutine ..." to get name list.

    if (/^\s*subroutine /i || /^\s*recursive *subroutine /i || 
        /^\s*function /i   || /^\s*recursive *function /i || 
        /^\s*integer *function /i || /^\s*real\(rp\) *function /i || 
        /^\s*logical *function /i || 
        /^\s*elemental *subroutine /i || /^\s*interface /i) {
      $name = $';              #' strip off "subroutine"
      $name =~ s/ *\(.*//;     # strip off " (..."
      $name =~ s/ +$//;        # strip off trailing blank if no arg list
      $name =~ tr/A-Z/a-z/; #lowercase
      chomp $name;
      # print "Hashed: $name\n";
     if (! exists $f90_hash{$name}) { $f90_hash{$name} = $file; }

      # skip rest of routine including contained routines

      $count = 1;
      while (<F_IN>) {
        if (/^\s*end /i) {
          $_ = $';  #'
          if (/^\s*subroutine/i || /^\s*function/i || /^\s*interface/i) {
            $count = $count - 1;
          }
        }
        elsif (/^\s*subroutine /i || /^\s*recursive subroutine /i || 
              /^\s*function /i || /^\s*elemental subroutine /i ||
              /^\s*real\(rp\) *function /i || /^\s*interface /i) {
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
