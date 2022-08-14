#!/usr/bin/perl

use File::Find;

#---------------------------------------------------------
# List of subroutines too low level to be mentioned

$pl_hash{"shortfactorial"} = "windowLS.f90";
$pl_hash{"fixedwindowls"} = "windowLS.f90";
$pl_hash{"initfixedwindowls"} = "windowLS.f90";
$pl_hash{"destfixedwindowls"} = "windowLS.f90";
$pl_hash{"re_allocate2_complex"} = "re_allocate_mod.f90";
$pl_hash{"palett"} = "pgsubs.f90";
$pl_hash{"attribute_free1"} = "attribute_mod.f90";
$pl_hash{"attribute_free2"} = "attribute_mod.f90";
$pl_hash{"attribute_free3"} = "attribute_mod.f90";
$pl_hash{"parser_add_branch"} = "bmad_parser_mod.f90";
$pl_hash{"parser_set_attribute"} = "bmad_parser_mod.f90";
$pl_hash{"bp_set_ran_status"} = "bmad_parser_mod.f90";
$pl_hash{"re_allocate_logical2d"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate_string2d"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate_integer2d"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate_real2d"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate2_real"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate2_integer"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate2_string"} = "re_allocate_mod.f90";
$pl_hash{"re_allocate2_logical"} = "re_allocate_mod.f90";
$pl_hash{"form_digested_bmad_file_name"} = "bmad_parser_mod.f90";
$pl_hash{"add_this_multipass"} = "bmad_parser_mod.f90";
$pl_hash{"settable_dep_var_bookkeeping"} = "bmad_parser_mod.f90";
$pl_hash{"compute_slave_coupler"} = "bookkeeper_mod.f90";
$pl_hash{"makeup_multipass_slave"} = "bookkeeper_mod.f90";
$pl_hash{"ran_uniform_vector"} = "random_mod.f90";
$pl_hash{"integration_timer"} = "integration_timer_mod.f90";
$pl_hash{"super_mrqcof"} = "super_mrqmin_mod.f90";
$pl_hash{"super_gaussj"} = "super_mrqmin_mod.f90";
$pl_hash{"init_attribute_name_array"} = "attribute_mod.f90";
$pl_hash{"parser_add_constant"} = "bmad_parser_mod.f90";
$pl_hash{"quad_mat2_calc"} = "make_mat6_mod.f90";
$pl_hash{"qrfac"} = "lmdif_mod.f90";
$pl_hash{"parser_add_lord"} = "bmad_parser_mod.f90";
$pl_hash{"mat4_multipole"} = "make_mat6_mod.f90";
$pl_hash{"evaluate_logical"} = "bmad_parser_mod.f90";
$pl_hash{"get_sequence_args"} = "bmad_parser_mod.f90";
$pl_hash{"allocate_plat"} = "bmad_parser_mod.f90";
$pl_hash{"add_taylor_term"} = "bmad_parser_mod.f90";
$pl_hash{"parser_file_stack"} = "bmad_parser_mod.f90";
$pl_hash{"get_next_word"} = "bmad_parser_mod.f90";
$pl_hash{"get_overlay_group_names"} = "bmad_parser_mod.f90";
$pl_hash{"make_this_overlay_group_lord"} = "bmad_parser_mod.f90";
$pl_hash{"parse_integer_list"} = "bmad_parser_mod.f90";
$pl_hash{"parser_read_old_format_sr_wake"} = "bmad_parser_mod.f90";
$pl_hash{"expect_one_of"} = "bmad_parser_mod.f90";
$pl_hash{"get_switch"} = "bmad_parser_mod.f90";
$pl_hash{"expect_this"} = "bmad_parser_mod.f90";
$pl_hash{"parse_real_list"} = "bmad_parser_mod.f90";
$pl_hash{"parser_identify_fork_to_element"} = "bmad_parser_mod.f90";
$pl_hash{"parser_get_logical"} = "bmad_parser_mod.f90";
$pl_hash{"parse_real_lists"} = "bmad_parser_mod.f90";
$pl_hash{"bmad_parser_string_attribute_set"} = "bmad_parser_mod.f90";
$pl_hash{"parser_error"} = "bmad_parser_mod.f90";
$pl_hash{"parser2_add_superimpose"} = "bmad_parser_mod.f90";
$pl_hash{"parser_init_custom_elements"} = "bmad_parser_mod.f90";
$pl_hash{"evaluate_array_index"} = "bmad_parser_mod.f90";
$pl_hash{"equal_sign_here"} = "bmad_parser_mod.f90";
$pl_hash{"parser_read_old_format_lr_wake"} = "bmad_parser_mod.f90";
$pl_hash{"check_for_superimpose_problem"} = "bmad_parser_mod.f90";
$pl_hash{"ix_far_index"} = "bmad_parser_mod.f90";
$pl_hash{"drift_multipass_name_correction"} = "bmad_parser_mod.f90";
$pl_hash{"parser_add_superimpose"} = "bmad_parser_mod.f90";
$pl_hash{"parser_read_lr_wake"} = "bmad_parser_mod.f90";
$pl_hash{"parse_evaluate_value"} = "bmad_parser_mod.f90";
$pl_hash{"parse_integer_list2"} = "bmad_parser_mod.f90";
$pl_hash{"parse_cylindrical_map"} = "bmad_parser_mod.f90";
$pl_hash{"parse_real_list2"} = "bmad_parser_mod.f90";
$pl_hash{"parse_taylor_field"} = "bmad_parser_mod.f90";
$pl_hash{"parser_get_integer"} = "bmad_parser_mod.f90";
$pl_hash{"get_list_of_names"} = "bmad_parser_mod.f90";
$pl_hash{"parse_grid_field"} = "bmad_parser_mod.f90";
$pl_hash{"parse_cartesian_map"} = "bmad_parser_mod.f90";
$pl_hash{"parser_read_sr_wake"} = "bmad_parser_mod.f90";
$pl_hash{"parser_print_line"} = "bmad_parser_mod.f90";
$pl_hash{"parse_line_or_list"} = "bmad_parser_mod.f90";
$pl_hash{"load_parse_line"} = "bmad_parser_mod.f90";
$pl_hash{"qromb_rad_int"} = "rad_int_common.f90";
$pl_hash{"init_bmad_parser_common"} = "bmad_parser_mod.f90";
$pl_hash{"bbi_slice_calc"} = "make_mat6_mod.f90";
$pl_hash{"reallocate_control"} = "bmad_utils_mod.f90";
$pl_hash{"compute_super_lord_s"} = "bmad_parser_mod.f90";
$pl_hash{"save_a_step"} = "em_field_mod.f90";
$pl_hash{"verify_valid_name"} = "bmad_parser_mod.f90";
$pl_hash{"get_called_file"} = "bmad_parser_mod.f90";
$pl_hash{"mexp"} = "bmad_basic_mod.f90";
$pl_hash{"integration_timer_fibre"} = "integration_timer_mod.f90";
$pl_hash{"propagate_part_way"} = "rad_int_common.f90";
$pl_hash{"reallocate_bp_com_const"} = "bmad_parser_mod.f90";
$pl_hash{"parser_debug_print_info"} = "bmad_parser_mod.f90";
$pl_hash{"parser_expand_line"} = "bmad_parser_mod.f90";
$pl_hash{"word_to_value"} = "bmad_parser_mod.f90";
$pl_hash{"term_diff"} = "integration_timer_mod.f90";
$pl_hash{"integration_timer_ele"} = "integration_timer_mod.f90";
$pl_hash{"word_read"} =               "word_read.f90";
$pl_hash{"calc_file_number"} =        "calc_file_number.f90";
$pl_hash{"fdjac2"} =                  "lmdif_mod.f90";
$pl_hash{"increment_file_number"} =   "increment_file_number.f90";
$pl_hash{"get_file_number"} =         "get_file_number.f90";
$pl_hash{"doubleup_quotes"} =         "doubleup_quotes.f90";
$pl_hash{"qrsolv"} =                  "lmdif_mod.f90";
$pl_hash{"change_file_number"} =      "change_file_number.f90";
$pl_hash{"file_get_open"} =           "file_get_open.f90";
$pl_hash{"file_directorizer"} =       "file_directorizer.f90";
$pl_hash{"make_legal_comment"} =      "make_legal_comment.f90";
$pl_hash{"ran_gauss_vector"} = "random_mod.f90";
$pl_hash{"probability_funct"} = "probability_funct.f90";
$pl_hash{"re_associate_real"} = "reallocate_mod.f90";
$pl_hash{"run_timer"} = "run_timer.f90";
$pl_hash{"str_upcase"} = "str_upcase.f90";
$pl_hash{"mat6_add_pitch"} = "bmad_utils_mod.f90";
$pl_hash{"re_associate_string"} = "reallocate_mod.f90";
$pl_hash{"sbend_body_with_k1_map"} = "make_mat6_mod.f90";
$pl_hash{"fullfilename"} = "fullfilename.f90";
$pl_hash{"re_allocate_real"} = "reallocate_mod.f90";
$pl_hash{"i_size"} = "bmad_and_cpp.f90";
$pl_hash{"lmpar"} = "lmdif_mod.f90";
$pl_hash{"ety"} = "eigen_mod.f90";
$pl_hash{"save_bunch_track"} = "save_bunch_track.f90";
$pl_hash{"value_to_line"} = "io_mod.f90";
$pl_hash{"eigensys"} = "eigen_mod.f90";
$pl_hash{"re_allocate_string"} = "reallocate_mod.f90";
$pl_hash{"query_string"} = "query_string.f90";
$pl_hash{"add_this_taylor_term"} = "bmad_parser_mod.f90";
$pl_hash{"r_size"} = "bmad_and_cpp.f90";
$pl_hash{"operator"} = "equality_mod.f90";
$pl_hash{"integrand"} = "ibs_mod.f90";
$pl_hash{"exp_bessi0"} = "touschek_mod.f90";
$pl_hash{"re_associate_logical"} = "reallocate_mod.f90";
$pl_hash{"eval_logical"} = "utilities_mod.f90";
$pl_hash{"file_get"} = "file_get.f90";
$pl_hash{"lmdif"} = "lmdif_mod.f90";
$pl_hash{"csr_bell"} = "csr_bell.f90";
$pl_hash{"word_len"} = "word_len.f90";
$pl_hash{"ety2"} = "eigen_mod.f90";
$pl_hash{"etyt"} = "eigen_mod.f90";
$pl_hash{"etdiv"} = "eigen_mod.f90";
$pl_hash{"inverse"} = "inverse.f90";
$pl_hash{"re_allocate_integer"} = "reallocate_mod.f90";
$pl_hash{"calc_wiggler_g_params"} = "rad_int_common.f90";
$pl_hash{"re_allocate_logical"} = "reallocate_mod.f90";
$pl_hash{"mat6_add_offsets"} = "mat6_add_offsets.f90";
$pl_hash{"assignment"} = "bmad_taylor_mod.f90";
$pl_hash{"re_associate_integer"} = "reallocate_mod.f90";

#---------------------------------------------------------
# make a list of names from bmad_subroutines.html

    #$name =~ s/protect\\parbox\{6in\}\{//;
    #$name =~ tr/A-Z/a-z/;    # lowercase
    #$name =~ s/\\//g;       # remove "\"
    #chomp $name;

$tex_file =  'subroutines.tex';
open (F_IN, $tex_file) || die ("Cannot open File: $tex_file");

while (<F_IN>) {
  if (/\\index\[routine\]\{(.*?)\}/i)  {      # match to "\index[routine]{...}" 
    $name = $1;
    if (exists $tex_hash{$name}) {next;}

    $_ = <F_IN>;
    if (!/\\item/) {$_ = <F_IN>;}
    if (!/\\item/) {
      print "No item for: $name\n";
      $tex_hash{$name} = $name;
      next;
    }

    while () {
      if (/\(/) {last;}  # Next line if no "("
      $_ = <F_IN>;
    }
    $arg_list = $_;
    $arg_list =~ s/.*?\(//;  # Strip off "...("
    chomp $arg_list;

    while () {
      if ($arg_list =~ /\)/) {   # If matches ")"
        $arg_list =~ s/\).*//;   # String ")..."
        last;
      }
      $_ = <F_IN>;
      $arg_list = $arg_list . $_;
      chomp $arg_list;
    }

    $arg_list =~ s/\\\\//g;           # Remove "\\"
    $arg_list =~ s/\\hspace.+?\}//g;  # Remove "\hspace{...}"
    $arg_list =~ s/\\hspace//g;       # Remove "\hfill"
    $arg_list =~ s/ {2,}/ /g;         # Multiple space to single space
    $tex_hash{$name} = $arg_list;
    ## print "$arg_list\n";
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
find (\&searchit, '../../sim_utils');

#---------------------------------------------------------
# find all BMAD subroutines and see if they are in the tex file

print "\n---------------------------------------------------\n";
print "Subroutines in bmad but not in subroutines.tex:\n\n";

while (($key, $val) = each %actual_arg_hash) {
  if (! exists $tex_hash{$key} && ! exists $pl_hash{$key} && ! exists $private_hash{$key}) {
#    printf "%-30s%s\n", ":$key:", $val;
    printf "\$pl_hash\{\"$key\"\} = \"$f90_hash{$key}\"\;\n";
  }
}

#---------------------------------------------------------
# find all listings in the tex file and see if there are bmad subroutines

print "\n---------------------------------------------------\n";
print "Subroutines in subroutines.tex that are not in bmad:\n\n";

while (($key, $val) = each %tex_hash) {
  if (! exists $actual_arg_hash{$key} || exists $private_hash{$key}) {
    print ":$key:\n";
  }
}

#---------------------------------------------------------
# find all listings in the tex file and see if there are bmad subroutines

print "\n---------------------------------------------------\n";
print "Subroutines in bmad_subroutines.pl that are not in bmad:\n\n";

while (($key, $val) = each %pl_hash) {
  if (! exists $actual_arg_hash{$key} || exists $private_hash{$key}) {
    print "\"$key\"\n";
  }
}

#---------------------------------------------------------
# Compare arg lists.

print "\n---------------------------------------------------\n";
print "Arg lists that are different in bmad and in subroutines.tex:\n\n";

while (($key, $actual_arg) = each %actual_arg_hash) {
  if (exists $private_hash{$key}) {next;}
  if (!exists $tex_hash{$key}) {next;}
  $tex_val = $tex_hash{$key};
  if (lc($tex_val) eq lc($actual_arg)) {next;}
  if ($actual_arg eq "") {next;}  # Ignore if f90 arg list is blank.
  printf "\n";
  printf "$key\n";
  printf "    f90: $actual_arg\n";
  printf "    tex: $tex_val\n";
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

    # record private routines

    if (/^ *private */i) {
      $this = $';   #' strip off "private"
      @names = split(/ *, */, $this);
      foreach $nam (@names) {
        chomp $nam;
        $private_hash{$nam} = "";
      }
      next;
    }

    # skip interface blocks

    if (/^ *interface *$/i) {    
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
        if (!exists $tex_hash{$name} && !exists $pl_hash{$name} && !exists $private_hash{$name}) {
          print "\nFile: $file\n";
          $this2 = $this; 
          $this2 =~ s/\s*\(.*//;
          print "\\index\[routine\]\{$this2\}\n";
          $this2 =~ s/_/./g;
          print "\\label\{r:$this2\}\n";
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

    if (/^\s*subroutine /i || /^\s*recursive *subroutine /i || /^\s*interface /i || 
        /^\s*function /i   || /^\s*recursive *function /i || 
        /^\s*integer *function /i || /^\s*real\(rp\) *function /i || 
        /^\s*logical *function /i || /^\s*pure *subroutine /i || 
        /^\s*elemental *subroutine /i || /^\s*elemental *function /i) {
      $name = $';              #' strip off "subroutine"
      $name =~ s/ *\(.*//;     # strip off " (..."
      $name =~ s/ +$//;        # strip off trailing blank if no arg list
      $name =~ tr/A-Z/a-z/; #lowercase
      chomp $name;

      if (/\((.*)/) {
        $arg_list = $1;
        while () {
          if ($arg_list =~ /\)/) {   # If matches ")"
            $arg_list =~ s/\).*//;   # String ")..."
            last;
          }
          $_ = <F_IN>;
          $arg_list = $arg_list . $_;
          chomp $arg_list;
        }
      } else {
        $arg_list = "";
      }

      $arg_list =~ s/\&//g;           # Remove "&"
      $arg_list =~ s/ {2,}/ /g;         # Multiple space to single space
        
      # print "Hashed: $name\n";
     if (! exists $actual_arg_hash{$name}) { 
       $actual_arg_hash{$name} = $arg_list; 
       $f90_hash{$name} = $file;
     }

      # skip rest of routine including contained routines

      $count = 1;
      while (<F_IN>) {
        if (/^\s*end /i) {
          $_ = $';  #'
          if (/^\s*subroutine/i || /^\s*function/i || /^\s*interface/i) {
            $count = $count - 1;
          }
        }
        elsif (/^\s*subroutine /i || /^\s*recursive *subroutine /i || 
              /^\s*function /i || /^\s*elemental *subroutine /i || /^\s*elemental *function /i || 
              /^\s*pure *subroutine /i || /^\s*real\(rp\) *function /i || /^\s*interface /i) {
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
