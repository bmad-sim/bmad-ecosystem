#!/usr/bin/perl

use File::Find;

#---------------------------------------------------------
# List of subroutines too low level to be mentioned

$tex_hash{"coord_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"amode_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"param_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"ran_uniform_vector"} = "random_mod.f90";
$tex_hash{"ran_uniform_scaler"} = "random_mod.f90";
$tex_hash{"integration_timer"} = "integration_timer_mod.f90";
$tex_hash{"read_xsif_wake"} = "xsif_parser.f90";
$tex_hash{"sr_wake_to_c"} = "bmad_and_cpp.f90";
$tex_hash{"test_f_em_field"} = "test_f_side.f90";
$tex_hash{"test_f_modes"} = "test_f_side.f90";
$tex_hash{"test_f_amode"} = "test_f_side.f90";
$tex_hash{"bmad_com_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"mode_info_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"lr_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"apply_wig_exp_int_ay"} = "symp_lie_mod.f90";
$tex_hash{"lr_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"test_f_floor_position"} = "test_f_side.f90";
$tex_hash{"sr_wake_in_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"coord_in_orbit_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"sr_wake_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"modes_to_f2"} = "bmad_and_cpp.f90";
$tex_hash{"test_f_ele"} = "test_f_side.f90";
$tex_hash{"test_f_wig_term"} = "test_f_side.f90";
$tex_hash{"test_f_control"} = "test_f_side.f90";
$tex_hash{"test_f_sr_wake"} = "test_f_side.f90";
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
$tex_hash{"gauss_int_liar"} = "macroparticle_mod.f90";
$tex_hash{"end_spline_calc"} = "spline_akima.f90";
$tex_hash{"track_it"} = "symp_lie_mod.f90";
$tex_hash{"parser_add_variable"} = "bmad_parser_mod.f90";
$tex_hash{"quad_mat2_calc"} = "make_mat6_mod.f90";
$tex_hash{"mat6_add_multipoles_and_s_offset"} = "make_mat6_bmad.f90";
$tex_hash{"qp_set_this_char_size"} = "quick_plot.f90";
$tex_hash{"simplify_path"} = "read_digested_bmad_file.f90";
$tex_hash{"write_this"} = "io_mod.f90";
$tex_hash{"write_out"} = "io_mod.f90";
$tex_hash{"str"} = "io_mod.f90";
$tex_hash{"rchomp"} = "io_mod.f90";
$tex_hash{"gauss_int"} = "macroparticle_mod.f90";
$tex_hash{"z_pitch_correction"} = "track1_bmad.f90";
$tex_hash{"qrfac"} = "lmdif_mod.f90";
$tex_hash{"slope_calc"} = "spline_akima.f90";
$tex_hash{"if_error"} = "if_error.f90";
$tex_hash{"track_bend_edge"} = "track1_mod.f90";
$tex_hash{"ion_kick_2d"} = "ion_kick_2d.f90";
$tex_hash{"twiss_to_taylor"} = "xsif_parser.f90";
$tex_hash{"add_ele"} = "xsif_parser.f90";
$tex_hash{"read_wake"} = "xsif_parser.f90";
$tex_hash{"add_t_term"} = "xsif_parser.f90";
$tex_hash{"xsif_error"} = "xsif_parser.f90";
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
$tex_hash{"from_action"} = "convert_coords.f90";
$tex_hash{"track_fwd"} = "mat627_mod.f90";
$tex_hash{"non_db_set"} = "bmad_to_db.f90";
$tex_hash{"allocate_saved_orbit"} = "em_field_mod.f90";
$tex_hash{"preparse_element_init"} = "bmad_parser_mod.f90";
$tex_hash{"da_z_dx__dx"} = "symp_lie_mod.f90";
$tex_hash{"da_z_dy"} = "symp_lie_mod.f90";
$tex_hash{"reallocate_control_"} = "bmad_utils_mod.f90";
$tex_hash{"da_z_dx__dy"} = "symp_lie_mod.f90";
$tex_hash{"compute_super_lord_s"} = "bmad_parser_mod.f90";
$tex_hash{"save_a_step"} = "em_field_mod.f90";
$tex_hash{"end_z_calc"} = "track1_bmad.f90";
$tex_hash{"verify_valid_name"} = "bmad_parser_mod.f90";
$tex_hash{"warning"} = "bmad_parser_mod.f90";
$tex_hash{"cache_quad_bend"} = "radiation_integrals.f90";
$tex_hash{"calc_g"} = "radiation_mod.f90";
$tex_hash{"zero_ave"} = "twiss_at_element.f90";
$tex_hash{"a_y__dx"} = "symp_lie_mod.f90";
$tex_hash{"a_y__dy"} = "symp_lie_mod.f90";
$tex_hash{"to_action"} = "convert_coords.f90";
$tex_hash{"delete_last_chars"} = "add_superimpose.f90";
$tex_hash{"a_y"} = "symp_lie_mod.f90";
$tex_hash{"bbi_kick_matrix"} = "make_mat6_bmad.f90";
$tex_hash{"b_ave"} = "quad_beta_ave.f90";
$tex_hash{"increment_pointer"} = "bmad_parser_mod.f90";
$tex_hash{"track_solenoid_edge"} = "boris_mod.f90";
$tex_hash{"get_called_file"} = "bmad_parser_mod.f90";
$tex_hash{"mexp"} = "bmad_basic_mod.f90";
$tex_hash{"integration_timer_fibre"} = "integration_timer_mod.f90";
$tex_hash{"dint_a_y_dx"} = "symp_lie_mod.f90";
$tex_hash{"db_init_it"} = "bmad_to_db.f90";
$tex_hash{"propagate_part_way"} = "rad_int_common.f90";
$tex_hash{"add_all_superimpose"} = "bmad_parser_mod.f90";
$tex_hash{"track_back"} = "mat627_mod.f90";
$tex_hash{"calc_g_params"} = "rad_int_common.f90";
$tex_hash{"diff"} = "integration_timer_mod.f90";
$tex_hash{"track1_sr_trans_wake"} = "macro_particle_mod.f90";
$tex_hash{"bmad_to_cesr_err_type"} = "bmad_to_cesr.f90";
$tex_hash{"multi_turn_func"} = "multi_turn_tracking_to_mat.f90";
$tex_hash{"evaluate_value"} = "bmad_parser_mod.f90";
$tex_hash{"dint_a_y_dx__dx"} = "symp_lie_mod.f90";
$tex_hash{"track_period"} = "track1_wiedemann_wiggler.f90";
$tex_hash{"dint_a_y_dx__dy"} = "symp_lie_mod.f90";
$tex_hash{"reallocate_bp_com_var"} = "bmad_parser_mod.f90";
$tex_hash{"error_exit"} = "bmad_parser_mod.f90";
$tex_hash{"delete_double_slash"} = "add_superimpose.f90";
$tex_hash{"do_vsp_eles"} = "create_vsp_volt_elements.f90";
$tex_hash{"word_to_value"} = "bmad_parser_mod.f90";
$tex_hash{"compute2_super_lord_s"} = "bmad_parser_mod.f90";
$tex_hash{"twiss_ave"} = "twiss_at_element.f90";
$tex_hash{"get_ele_theory"} = "bmad_to_db.f90";
$tex_hash{"makeup_overlay_slave"} = "bookkeeper_mod.f90";
$tex_hash{"da_z_dy__dx"} = "symp_lie_mod.f90";
$tex_hash{"da_z_dy__dy"} = "symp_lie_mod.f90";
$tex_hash{"rkck_bmad"} = "runge_kutta_mod.f90";
$tex_hash{"bookit"} = "create_group.f90";
$tex_hash{"term_diff"} = "integration_timer_mod.f90";
$tex_hash{"track1_bunch"} = "macro_particle_mod.f90";
$tex_hash{"insert_info"} = "bmad_to_cesr.f90";
$tex_hash{"pushit"} = "bmad_parser_mod.f90";
$tex_hash{"type_get"} = "bmad_parser_mod.f90";
$tex_hash{"da_z_dx"} = "symp_lie_mod.f90";
$tex_hash{"get_taylor"} = "integration_timer_mod.f90";
$tex_hash{"sr_long_wake_calc"} = "macro_particle_mod.f90";
$tex_hash{"order_macroparticles_in_z"} = "macro_particle_mod.f90";
$tex_hash{"integration_timer_ele"} = "integration_timer_mod.f90";
$tex_hash{"map_index"} = "ptc_interface_mod.f90";
$tex_hash{"solenoid_mat_calc"} = "make_mat6_mod.f90";

$tex_hash{"coef23_calc"} =             "spline_akima.f90";
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
$tex_hash{"load_data"} =               "quick_plot.f90";
$tex_hash{"spawn_command"} =           "spawn_command.f90";
$tex_hash{"to_node_name"} =            "to_node_name.f90";
$tex_hash{"file_get_open"} =           "file_get_open.f90";
$tex_hash{"match_word_exact"} =        "match_word_exact.f90";
$tex_hash{"file_directorizer"} =       "file_directorizer.f90";
$tex_hash{"make_legal_comment"} =      "make_legal_comment.f90";
$tex_hash{"ran_gauss_vector"} = "random_mod.f90";
$tex_hash{"ran_gauss_scaler"} = "random_mod.f90";
$tex_hash{"plot_it"} = "plot_example.f90";
$tex_hash{"update_wig_coefs"} = "symp_lie_mod.f90";
$tex_hash{"compute_slave_aperture"} = "bookkeeper_mod.f90";
$tex_hash{"update_wig_y_terms"} = "symp_lie_mod.f90";
$tex_hash{"update_wig_x_s_terms"} = "symp_lie_mod.f90";
$tex_hash{"bsq_drift1"} = "symp_lie_mod.f90";
$tex_hash{"bsq_drift2"} = "symp_lie_mod.f90";
$tex_hash{"zero_this_track"} = "track_many.f90";
$tex_hash{"apply_p_x"} = "symp_lie_mod.f90";
$tex_hash{"apply_p_y"} = "symp_lie_mod.f90";


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
    }

    $old = $_;
  }

  close (F_IN);

  $_ = $file;

}
