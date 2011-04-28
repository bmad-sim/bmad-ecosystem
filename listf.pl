#!/usr/bin/perl
# modify to use File::Spec for transparent use on VMS
# 2006.01.18 mjf

use File::Find;
use File::Spec::Functions;

# Following two lines should be uncommented on VMS to check the master_list:

#my @args = ("list_list", @ARGV);
#system(@args) == 0 or warn "system @args failed: $!\n";

$found_one = 0;

my $curdir = curdir();
my $updir  = updir();

if (-r catfile( $curdir, "bmad", "modules", "bmad_struct.f90" )) {
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

if (-r catfile( $curdir, "cesr_utils", "modules", "cesr_utils.f90" )) {
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


if (-r catfile( $curdir, "sim_utils", "interfaces", "sim_utils.f90" )) {
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

if (-r catfile( $curdir, "mpm_utils", "code", "butout.f90" )) {
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


if (-r catfile( $curdir, "recipes_f-90_LEPP", "lib_src", "nr.f90" )) {
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


if (-r catfile( $curdir, "forest", "code", "i_tpsa.f90" )) {
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

if (-r catfile( $curdir, "tao", "code", "tao_struct.f90" )) {
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

if (-r catfile( $curdir, "bmadz", "modules", "bmadz_struct.f90" )) {
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

if (-r catfile( $curdir, "nonlin_bpm", "code", "nonlin_bpm_init.f90" )) {
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

# Look for arguments

$extra = 0;
$also_match_params = 1;

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

  last;

}

#

$match_str = @field[$i];
$match_str =~ s/\*/\\w\*/g;  # replace "*" by "\w*"

$match_param_str = $match_str;
if ($match_param_str =~ /\\w\*$/) {$match_param_str = $match_param_str . '\\$*';}
if ($match_param_str =~ /\$$/) {$match_param_str =~ s/\$$/\\\$/;}

# Help if needed

if ($match_str eq "") {
  print "Usage:\n";
  print "  listf {options} <search_string>\n\n";
  print "Options:\n";
  print "   -d <dir>    # Search files in <dir> and sub-directories for matches.\n";
  ## print "   -p          # Match to parameters as well.\n";
  exit;
}

# Search for a match

find(\&searchit, $bmad_dir);
find(\&searchit, $sim_utils_dir);
find(\&searchit, $cesr_utils_dir);
find(\&searchit, $tao_dir);
find(\&searchit, $mpm_utils_dir);
find(\&searchit, $bmadz_dir);
find(\&searchit, $nonlin_bpm_dir);
if ($extra == 1) {find(\&searchit, $extra_dir);}
## find(\&searchit, $recipes_dir);
## find(\&searchit, $forest_dir);

if ($found_one == 0) {
  print "Cannot match string! $match_str\n";
  ## if ($also_match_params == 1) {print "Also cannot match parameter string $match_param_str\n";}
  print "Note: getf does not search numerical recipes, forest, nor any program directories.\n";
  ## if ($also_match_params == 0) {print "Note: Use '-p' option to search for parameters as well as routines and structures.\n";}
} else {
  print "\n";
}

#---------------------------------------------------------

sub searchit {

  if ($File::Find::name =~ /cesr_utils\/f77/) {return;}  # Do not search this directory.

  if (/\#/) {return;}
  $file = $_;

# If a f90 file look in the file for a match.

  if ($file =~ /\.f90$/i) {

    $found_in_file = 0;
    $in_module_header = 0;

    $stat = open (F_IN, $file); 
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    while (<F_IN>) {

      if (/^ *interface *$/i) {   # skip interface blocks
        while (<F_IN>) {
          if (/^ *end +interface/i) {last;}
        }
      }

      if (/^ *type *\(/i) {next;}   # skip "type (" constructs

      # is this a module?

      if (/^ *module /i) {
        $in_module_header = 1;
      }

      # is this a parameter?

      if ($also_match_params ==1 && $in_module_header ==1 && /, *parameter *::/i) {
        $rest_of_line = $';              #' grab the parameter name
        $line = $_;
        $param_match = 0;
        while (1) {
          ## print "1: $rest_of_line\n";
          if ($rest_of_line eq "") {last;}
          if (! ($rest_of_line =~ /^\s*(\w+\$*) *([\(=])/)) {last;} 
          $param_name = $1;
          $rest_of_line = $';      #'
          if ($2 eq "\(") {$rest_of_line =~ s/^.*?\) *=//;}
          ## print "2: $param_name &&& $rest_of_line\n";
          if ($param_name =~ /^$match_param_str$/i) {
            $param_match = 1;
            last;
          }
          # strip off param value
          if ($rest_of_line =~ /^ *\(/) {          # For "param = (/ ... /)," construct
            $rest_of_line =~ s/^ *\(.*?\) *, *//; 
          } elsif ($rest_of_line =~ /^ *\[/) {     # For "param = [ ... ]," construct     
            $rest_of_line =~ s/^ *\[.*?\] *, *//;
          } else {
            $rest_of_line =~ s/^.*?, *//;
          }

        }

        if ($param_match == 1) {
          if ($file eq "$_\.f90") {$matches_file_name = 1;}
          $found_one = 1;
          if ($found_in_file == 0) {print "File: $File::Find::name\n";}
          print "  $line";
        }
      }

      # if a routine then look for match

      if (&routine_here) {

        # If a match then print info

        if ($routine_name =~ /^\s*$match_str[ \(\n]/i) {
          if ($found_in_file == 0) {print "File: $File::Find::name\n";}
          $found_one = 1;
          $found_in_file = 1;
          print "     $_";
          if (/&\s*$/) {
            $_ = <F_IN>;
            print "   $_";
          }
        }

        # skip rest of routine including contained routines

        $count = 1;
        while (<F_IN>) {
          if (/^ *type *\(/i) {next;}   # ignore "type (" constructs
          if (/^ *end/i) {
            $_ = $';  #'
            if (/^ *subroutine/i || /^ *function/i ||
                /^ *type/i || /^ *interface/i) {$count = $count - 1;}
          }
          elsif (&routine_here) {
            $count = $count + 1;
          }
          if ($count == 0) {last;}
        }

      }
    }
    close (F_IN);

# match to C++ files

  } elsif ($file =~ /\.cpp$/ || $file =~ /\.h$/) {


    $found_in_file = 0;

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    while (<F_IN>) {

      if (/^ *class +(\w+) +\{ *$/) {     # match to "class xyz {"
        if ($1 =~ /^$match_str$/i) {
          $found_one = 1;
          if (!$found_in_file) {print "\n$File::Find::name\n";}

          $found_in_file = 1;
          print $_;
        }
      }
    }
    close (F_IN);
  }
  
# See if the file name matches.

  if (/^$match_str\.f90$/ && $found_in_file == 0) {
    $found_one = 1;
    print "$File::Find::name\n";
  }

#

  $_ = $file;

}

#---------------------------------------------------------

sub routine_here {

  if (/^\s*subroutine /i || /^\s*recursive subroutine /i || 
      /^\s*elemental subroutine /i ||
      /^\s*function /i || /^\s*recursive function /i ||
      /^\s*real\(rp\) *function /i || /^\s*integer *function /i ||
      /^\s*logical *function /i || /^\s*type /i || /^\s*interface /i) {
    $routine_name = $';              #' strip off "routine" string
    return 1;
  }

  return 0;

}
