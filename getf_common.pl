sub setup_dirs {

  # The idea is to look for a local copy of the library to search.
  # We have found a local copy when we find one specific file that we know 
  # is in the library.
  # -r is a file test operator which checks if the file is can be read

  if (-r catfile( $curdir, "bmad", "modules", "bmad_struct.f90")) {
    $bmad_dir = catfile( $curdir, "bmad" );
  } elsif (-r catfile( $updir, "bmad", "modules", "bmad_struct.f90")) {
    $bmad_dir = catfile( $updir, "bmad" );
  } elsif (-r catfile( $updir, $updir, "bmad", "modules", "bmad_struct.f90")) {
    $bmad_dir = catfile( $updir, $updir, "bmad" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "bmad", "modules", "bmad_struct.f90")) {
    $bmad_dir = catfile( $ENV{"ACC_SRC"}, "bmad" );
  } else {
    $bmad_dir = catfile( $ENV{"DIST_BASE_DIR"}, "bmad" );
  }

  if (-r catfile( $curdir, "cesr_utils", "modules", "cesr_utils.f90")) {
    $cesr_utils_dir = catfile( $curdir, "cesr_utils" );
  } elsif (-r catfile( $updir, "cesr_utils", "modules", "cesr_utils.f90")) {
    $cesr_utils_dir = catfile( $updir, "cesr_utils" );
  } elsif (-r catfile( $updir, $updir, "cesr_utils", "modules", "cesr_utils.f90")) {
    $cesr_utils_dir = catfile( $updir, $updir, "cesr_utils" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "cesr_utils")) {
    $cesr_utils_dir = catfile( $ENV{"ACC_SRC"}, "cesr_utils" );
  } else {
    $cesr_utils_dir = catfile( $ENV{"DIST_BASE_DIR"}, "cesr_utils" );
  }


  if (-r catfile( $curdir, "sim_utils", "interfaces", "sim_utils.f90")) {
    $sim_utils_dir = catfile( $curdir, "sim_utils" );
  } elsif (-r catfile( $updir, "sim_utils", "interfaces", "sim_utils.f90")) {
    $sim_utils_dir = catfile( $updir, "sim_utils" );
  } elsif (-r catfile( $updir, $updir, "sim_utils", "interfaces", "sim_utils.f90")) {
    $sim_utils_dir = catfile( $updir, $updir, "sim_utils" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "sim_utils")) {
    $sim_utils_dir = catfile( $ENV{"ACC_SRC"}, "sim_utils" );
  } else {
    $sim_utils_dir = catfile( $ENV{"DIST_BASE_DIR"}, "sim_utils" );
  }

  if (-r catfile( $curdir, "mpm_utils", "code", "butout.f90")) {
    $mpm_utils_dir = catfile( $curdir, "mpm_utils" );
  } elsif (-r catfile( $updir, "mpm_utils", "code", "butout.f90")) {
    $mpm_utils_dir = catfile( $updir, "mpm_utils" );
  } elsif (-r catfile( $updir, $updir, "mpm_utils", "code", "butout.f90")) {
    $mpm_utils_dir = catfile( $updir, $updir, "mpm_utils" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "mpm_utils")) {
    $mpm_utils_dir = catfile( $ENV{"ACC_SRC"}, "mpm_utils" );
  } else {
    $mpm_utils_dir = catfile( $ENV{"DIST_BASE_DIR"}, "mpm_utils" );
  }


  if (-r catfile( $curdir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
    $recipes_dir = catfile( $curdir, "recipes_f-90_LEPP" );
  } elsif (-r catfile( $updir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
    $recipes_dir = catfile( $updir, "recipes_f-90_LEPP" );
  } elsif (-r catfile( $updir, $updir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
    $recipes_dir = catfile( $updir, $updir, "recipes_f-90_LEPP" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "recipes_f-90_LEPP")) {
    $recipes_dir = catfile( $ENV{"ACC_SRC"}, "recipes_f-90_LEPP" );
  } else {
    $recipes_dir = catfile( $ENV{"DIST_BASE_DIR"}, "recipes_f-90_LEPP" );
  }


  if (-r catfile( $curdir, "forest", "code", "i_tpsa.f90")) {
    $forest_dir = catfile( $curdir, "forest" );
  } elsif (-r catfile( $updir, "forest", "code", "i_tpsa.f90")) {
    $forest_dir = catfile( $updir, "forest" );
  } elsif (-r catfile( $updir, $updir, "forest", "code", "i_tpsa.f90")) {
    $forest_dir = catfile( $updir, $updir, "forest" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "forest")) {
    $forest_dir = catfile( $ENV{"ACC_SRC"}, "forest" );
  } else {
    $forest_dir = catfile( $ENV{"DIST_BASE_DIR"}, "forest" );
  }

  if (-r catfile( $curdir, "tao", "code", "tao_struct.f90")) {
    $tao_dir = catfile( $curdir, "tao" );
  } elsif (-r catfile( $updir, "tao", "code", "tao_struct.f90")) {
    $tao_dir = catfile( $updir, "tao" );
  } elsif (-r catfile( $updir, $updir, "tao", "code", "tao_struct.f90")) {
    $tao_dir = catfile( $updir, $updir, "tao" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "tao")) {
    $tao_dir = catfile( $ENV{"ACC_SRC"}, "tao" );
  } else {
    $tao_dir = catfile( $ENV{"DIST_BASE_DIR"}, "tao" );
  }

  if (-r catfile( $curdir, "bmadz", "modules", "bmadz_struct.f90")) {
    $bmadz_dir = catfile( $curdir, "bmadz" );
  } elsif (-r catfile( $updir, "bmadz", "modules", "bmadz_struct.f90")) {
    $bmadz_dir = catfile( $updir, "bmadz" );
  } elsif (-r catfile( $updir, $updir, "bmadz", "modules", "bmadz_struct.f90")) {
    $bmadz_dir = catfile( $updir, $updir, "bmadz" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "bmadz")) {
    $bmadz_dir = catfile( $ENV{"ACC_SRC"}, "bmadz" );
  } else {
    $bmadz_dir = catfile( $ENV{"DIST_BASE_DIR"}, "bmadz" );
  }

  if (-r catfile( $curdir, "nonlin_bpm", "code", "nonlin_bpm_init.f90")) {
    $nonlin_bpm_dir = catfile( $curdir, "nonlin_bpm" );
  } elsif (-r catfile( $updir, "nonlin_bpm", "code", "nonlin_bpm_init.f90")) {
    $nonlin_bpm_dir = catfile( $updir, "nonlin_bpm" );
  } elsif (-r catfile( $updir, $updir, "nonlin_bpm", "code", "nonlin_bpm_init.f90")) {
    $nonlin_bpm_dir = catfile( $updir, $updir, "nonlin_bpm" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "nonlin_bpm")) {
    $nonlin_bpm_dir = catfile( $ENV{"ACC_SRC"}, "nonlin_bpm" );
  } else {
    $nonlin_bpm_dir = catfile( $ENV{"DIST_BASE_DIR"}, "nonlin_bpm" );
  }

  if (-r catfile( $curdir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
    $recipes_dir = catfile( $curdir, "recipes_f-90_LEPP" );
  } elsif (-r catfile( $updir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
    $recipes_dir = catfile( $updir, "recipes_f-90_LEPP" );
  } elsif (-r catfile( $updir, $updir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
    $recipes_dir = catfile( $updir, $updir, "recipes_f-90_LEPP" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "recipes_f-90_LEPP")) {
    $recipes_dir = catfile( $ENV{"ACC_SRC"}, "recipes_f-90_LEPP" );
  } else {
    $recipes_dir = catfile( $ENV{"DIST_BASE_DIR"}, "recipes_f-90_LEPP" );
  }

  if (-r catfile( $curdir, "forest", "code", "a_scratch_size.f90")) {
    $forest_dir = catfile( $curdir, "forest" );
  } elsif (-r catfile( $updir, "forest", "code", "a_scratch_size.f90")) {
    $forest_dir = catfile( $updir, "forest" );
  } elsif (-r catfile( $updir, $updir, "forest", "code", "a_scratch_size.f90")) {
    $forest_dir = catfile( $updir, $updir, "forest" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "forest")) {
    $forest_dir = catfile( $ENV{"ACC_SRC"}, "forest" );
  } else {
    $forest_dir = catfile( $ENV{"DIST_BASE_DIR"}, "forest" );
  }

  if (-r catfile( $curdir, "bsim", "code")) {
    $bsim_dir = catfile( $curdir, "bsim" );
  } elsif (-r catfile( $updir, "bsim", "code")) {
    $bsim_dir = catfile( $updir, "bsim" );
  } elsif (-r catfile( $updir, $updir, "bsim", "code")) {
    $bsim_dir = catfile( $updir, $updir, "bsim" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "bsim")) {
    $bsim_dir = catfile( $ENV{"ACC_SRC"}, "bsim" );
  } else {
    $bsim_dir = catfile( $ENV{"DIST_BASE_DIR"}, "bsim" );
  }

  if (-r catfile( $curdir, "bsim_cesr", "code")) {
    $bsim_cesr_dir = catfile( $curdir, "bsim_cesr" );
  } elsif (-r catfile( $updir, "bsim_cesr", "code")) {
    $bsim_cesr_dir = catfile( $updir, "bsim_cesr" );
  } elsif (-r catfile( $updir, $updir, "bsim_cesr", "code")) {
    $bsim_cesr_dir = catfile( $updir, $updir, "bsim_cesr" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "bsim_cesr")) {
    $bsim_cesr_dir = catfile( $ENV{"ACC_SRC"}, "bsim_cesr" );
  } else {
    $bsim_cesr_dir = catfile( $ENV{"DIST_BASE_DIR"}, "bsim_cesr" );
  }

  if (-r catfile( $curdir, "cesr_programs", "code")) {
    $cesr_programs_dir = catfile( $curdir, "cesr_programs" );
  } elsif (-r catfile( $updir, "cesr_programs", "code")) {
    $cesr_programs_dir = catfile( $updir, "cesr_programs" );
  } elsif (-r catfile( $updir, $updir, "cesr_programs", "code")) {
    $cesr_programs_dir = catfile( $updir, $updir, "cesr_programs" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "cesr_programs")) {
    $cesr_programs_dir = catfile( $ENV{"ACC_SRC"}, "cesr_programs" );
  } else {
    $cesr_programs_dir = catfile( $ENV{"DIST_BASE_DIR"}, "cesr_programs" );
  }

  if (-r catfile( $curdir, "cesrv", "code")) {
    $cesrv_dir = catfile( $curdir, "cesrv" );
  } elsif (-r catfile( $updir, "cesrv", "code")) {
    $cesrv_dir = catfile( $updir, "cesrv" );
  } elsif (-r catfile( $updir, $updir, "cesrv", "code")) {
    $cesrv_dir = catfile( $updir, $updir, "cesrv" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "cesrv")) {
    $cesrv_dir = catfile( $ENV{"ACC_SRC"}, "cesrv" );
  } else {
    $cesrv_dir = catfile( $ENV{"DIST_BASE_DIR"}, "cesrv" );
  }

  if (-r catfile( $curdir, "util_programs", "code")) {
    $util_programs_dir = catfile( $curdir, "util_programs" );
  } elsif (-r catfile( $updir, "util_programs", "code")) {
    $util_programs_dir = catfile( $updir, "util_programs" );
  } elsif (-r catfile( $updir, $updir, "util_programs", "code")) {
    $util_programs_dir = catfile( $updir, $updir, "util_programs" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "util_programs")) {
    $util_programs_dir = catfile( $ENV{"ACC_SRC"}, "util_programs" );
  } else {
    $util_programs_dir = catfile( $ENV{"DIST_BASE_DIR"}, "util_programs" );
  }

  if (-r catfile( $curdir, "examples", "simple_program")) {
    $examples_dir = catfile( $curdir, "examples" );
  } elsif (-r catfile( $updir, "examples", "simple_program")) {
    $examples_dir = catfile( $updir, "examples" );
  } elsif (-r catfile( $updir, $updir, "examples", "simple_program")) {
    $examples_dir = catfile( $updir, $updir, "examples" );
  } elsif (-r catfile( $ENV{"ACC_SRC"}, "examples")) {
    $examples_dir = catfile( $ENV{"ACC_SRC"}, "examples" );
  } else {
    $examples_dir = catfile( $ENV{"DIST_BASE_DIR"}, "examples" );
  }

  # Look for arguments

  $extra = 0;
  $also_match_params = 1; 
  $search_all = 0;

  $_ = "@ARGV";  
  @field = split;

  for ($i = 0; $i <= @field; $i++) {

    $_ = $field[$i];

    if (! /^\-/) {last;}

    if ($_ eq '-d') {
      $extra = 1; 
      $extra_dir = $field[$i+1];
      $i = $i + 1;
      next;
    }

    if ($_ eq '-p') {
      $also_match_params = 1; 
      next;
    }

    if ($_ eq '-a') {
      $search_all = 1;
      next;
    }

    last;

  }

}

#--------------------------------------------------------------------

sub search_it {

  my $curdir = curdir();
  my $updir = updir();

  #
   
  $match_str = $field[$i];
  $match_str =~ s/\*/\\w\*/g;   # replace "*" by "\w*" (any word characters)

  $match_param_str = $match_str;
  if ($match_param_str =~ /\\w\*$/) {$match_param_str = $match_param_str . '\\$*';}
  if ($match_param_str =~ /\$$/) {$match_param_str =~ s/\$$/\\\$/;}

  # Help if needed

  if ($match_str eq "") {
    print "Usage:\n";
    print "  @_ {options} <search_string>\n\n";
    print "Options:\n";
    print "   -a          # Search Numerical recipes, forest, and varies program directories as well.\n";
    print "   -d <dir>    # Search files in <dir> and sub-directories for matches.\n";
    ## print "   -p          # Match to parameters as well.\n";   # Not yet implemented!
    exit;
  }

  # Search for a match

  if ($extra == 1) {find(\&searchit, $extra_dir);}

  if ($search_all ==1) {
    find(\&searchit, $recipes_dir);
    find(\&searchit, $forest_dir);
    find(\&searchit, $bsim_dir);
    find(\&searchit, $bsim_cesr_dir);
    find(\&searchit, $cesr_programs_dir);
    find(\&searchit, $cesrv_dir);
    find(\&searchit, $util_programs_dir);
    find(\&searchit, $examples_dir);
  }

  find(\&searchit, $bmad_dir);
  find(\&searchit, $sim_utils_dir);
  find(\&searchit, $cesr_utils_dir);
  find(\&searchit, $tao_dir);
  find(\&searchit, $mpm_utils_dir);
  find(\&searchit, $bmadz_dir);
  find(\&searchit, $nonlin_bpm_dir);

  if ($found_one == 0) {
    print "Cannot match string! $match_str\n";
    ## if ($also_match_params == 1) {print "Also cannot match parameter string $match_param_str\n";}
    if ($search_all == 0) {
      print "Note: numerical recipes, forest, nor any program directories are searched.\n";
      print "      Use the '-a' option to expand the search.\n";
    }
    print "Note: C/C++ files are not searched.\n";
    if ($also_match_params == 0) {print "Note: Use '-p' option to search for parameters as well as routines and structures.\n";}
  } else {
    print "\n";
  }

}

#---------------------------------------------------------

sub routine_here {

  if (/^\s*subroutine /i || /^\s*recursive subroutine /i || 
      /^\s*elemental subroutine /i || /^\s*program /i || 
      /^\s*function /i || /^\s*recursive function /i ||
      /^\s*real\(rp\) *function /i || /^\s*integer *function /i ||
      /^\s*logical *function /i || /^\s*interface /i ||
      (/^\s*type /i && $_[0] == 1)) {
    $routine_name = $';              #' strip off "routine" string
    return 1;
  }

  return 0;

}

#-------------------------------------

$return = 1;   # So perl will not complain about not returning a true value

