#!/usr/bin/perl
# modify to use File::Spec for transparent use on VMS
# 2005.09.28 mjf

use File::Find;
use File::Spec::Functions;

# Following two lines should be uncommented on VMS to check the master_list:

#my @args = ("type_header", @ARGV);
#system(@args) == 0 or warn "system @args failed: $!\n";

$found_one = 0;

# The idea is to look for a local copy of the library to search.
# We have found a local copy when we find one specific file that we know 
# is in the library.
# -r is a file test operator which checks if the file is can be read

my $curdir = curdir();
my $updir = updir();

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
 
$match_str = $field[$i];
$match_str =~ s/\*/\\w\*/g;   # replace "*" by "\w*" (any word characters)

$match_param_str = $match_str;
if ($match_param_str =~ /\\w\*$/) {$match_param_str = $match_param_str . '\\$*';}
if ($match_param_str =~ /\$$/) {$match_param_str =~ s/\$$/\\\$/;}

# Help if needed

if ($match_str eq "") {
  print "Usage:\n";
  print "  getf {options} <search_string>\n\n";
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
  if ($also_match_params == 0) {print "Note: Use '-p' option to search for parameters as well as routines and structures.\n";}
} else {
  print "\n";
}

#---------------------------------------------------------

sub searchit {

  # Do not search this directory.
  if ($File::Find::name =~ /cesr_utils\/original_from_vms/) {return;}  

  if (/\#/) {return;}
  $file = $_;
  $matches_file_name = 0;

  # Check contents of .f90 file.

  if ($file =~ /\.f90$/) {
    $n_com = 0;
    $in_module_header = 0;

    @comments = ();     

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    while (<F_IN>) {

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
          print "\n$File::Find::name\n";
          print "$line";
        }
      }

      # skip interface blocks

      if (/^ *interface *$/i) {
        while (<F_IN>) {
          if (/^ *end +interface/i) {last;}
        }
      }

      # skip "type (" constructs and separator comments.

      if (/^ *type *\(/i) {next;}   
      if ($n_com == 0 && /^!\-\-\-\-\-\-\-\-\-/) {next;}
      if (/^#/) {next;}

      # Add to comment block if a comment line

      if (/^!/) {
        @comments[$n_com] = $_;
        $n_com++

      # match to type statement

      } elsif (/^ *type +$match_str[ \n]/i) {
        $found_one = 1;
        print "\n$File::Find::name\n";
        for ($i = 0; $i < $n_com; $i++) {print @comments[$i];}
        print "\n";
        print $_;
        while (<F_IN>) {
          print $_;
          if (/^ *end type/i) {last;}
        }

        $n_com = 0;

      # match to subroutine, function, etc.

      } elsif (&routine_here) {
        $in_module_header = 0;
        $routine_name =~ /^\s*(\w+)[ |\(\n]/i;
        $_ = $1;     # strip off "subroutine"
        if ($file eq "$_\.f90") {$matches_file_name = 1;}

        if (/^$match_str$/i) {
          $found_one = 1;
          print "\n$File::Find::name\n";
          for ($i = 0; $i < $n_com; $i++) {print @comments[$i];}
        }

        $n_com = 0;

        # skip rest of routine including contained routines

        $count = 1;
        while (<F_IN>) {
          if (/^ *end /i) {
            $_ = $';     #' 
            if (/^ *subroutine/i || /^ *function/i || /^ *interface/i) {
              $count = $count - 1;
            }
          }
          elsif (&routine_here) {
            $count = $count + 1;
          }
          if ($count == 0) {last;}
        }

      # If not a blank line, reset comment block

      } elsif (/[^\s]/) {
        $n_com = 0; 

      }

    }
    close (F_IN);
    
  # match to C++ files

  } elsif ($file =~ /\.cpp$/ || $file =~ /\.h$/) {

    $in_class = 0;

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    while (<F_IN>) {

      if (/^ *class +(\w+) +\{ *$/) {     # match to "class xyz {"
        if ($1 =~ /^$match_str$/i) {
          $in_class = 1;
          $found_one = 1;
          print "\n$File::Find::name\n";
        }
      }

      if ($in_class) {print $_;}
      if (m@// +end class@i) {$in_class = 0;}

    }
    close (F_IN);
  }

  #--------------------------------------------------
  # See if the file name matches.
  # Only print the file name here if:
  #    1) There are some comments at the top of the file or
  #    2) There has not been a match to a routine within the file.

  if ($matches_file_name == 0 && $file =~ /^$match_str\.f90$/i) {

    $found_one = 1;
    print "\n$File::Find::name\n";

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    # Print header

    $have_printed = 0;

    while (<F_IN>) {
      if (&routine_here) {last;}
      if (/^!/) {
        $have_printed = 1;
        print "$_";
      } elsif ($have_printed == 1) {
        last;
      }
    }

    close (F_IN);

  }

  $_ = $file;

}

#---------------------------------------------------------

sub routine_here {

  if (/^\s*subroutine /i || /^\s*recursive subroutine /i || 
      /^\s*elemental subroutine /i ||
      /^\s*function /i || /^\s*recursive function /i ||
      /^\s*real\(rp\) *function /i || /^\s*integer *function /i ||
      /^\s*logical *function /i || /^\s*interface /i) {
    $routine_name = $';              #' strip off "routine" string
    return 1;
  }

  return 0;

}
