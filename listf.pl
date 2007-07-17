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
} else {
  $bmad_dir = catfile( $ENV{"CESR_SRC"}, "bmad" );
}

if (-r catfile( $curdir, "cesr_utils", "modules", "cesr_utils.f90" )) {
  $cesr_utils_dir = catfile( $curdir, "cesr_utils" );
} elsif (-r catfile( $updir, "cesr_utils", "modules", "cesr_utils.f90")) {
  $cesr_utils_dir = catfile( $updir, "cesr_utils" );
} elsif (-r catfile( $updir, $updir, "cesr_utils", "modules", "cesr_utils.f90")) {
  $cesr_utils_dir = catfile( $updir, $updir, "cesr_utils" );
} else {
  $cesr_utils_dir = catfile( $ENV{"CESR_SRC"}, "cesr_utils" );
}


if (-r catfile( $curdir, "dcslib", "modules", "dcslib.f90" )) {
  $dcslib_dir = catfile( $curdir, "dcslib" );
} elsif (-r catfile( $updir, "dcslib", "modules", "dcslib.f90")) {
  $dcslib_dir = catfile( $updir, "dcslib" );
} elsif (-r catfile( $updir, $updir, "dcslib", "modules", "dcslib.f90")) {
  $dcslib_dir = catfile( $updir, $updir, "dcslib" );
} else {
  $dcslib_dir = catfile( $ENV{"CESR_SRC"}, "dcslib" );
}


if (-r catfile( $curdir, "recipes_f-90_LEPP", "lib_src", "nr.f90" )) {
  $recipes_dir = catfile( $curdir, "recipes_f-90_LEPP" );
} elsif (-r catfile( $updir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
  $recipes_dir = catfile( $updir, "recipes_f-90_LEPP" );
} elsif (-r catfile( $updir, $updir, "recipes_f-90_LEPP", "lib_src", "nr.f90")) {
  $recipes_dir = catfile( $updir, $updir, "recipes_f-90_LEPP" );
} else {
  $recipes_dir = catfile( $ENV{"CESR_SRC"}, "recipes_f-90_LEPP" );
}


if (-r catfile( $curdir, "forest", "code", "i_tpsa.f90" )) {
  $forest_dir = catfile( $curdir, "forest" );
} elsif (-r catfile( $updir, "forest", "code", "i_tpsa.f90")) {
  $forest_dir = catfile( $updir, "forest" );
} elsif (-r catfile( $updir, $updir, "forest", "code", "i_tpsa.f90")) {
  $forest_dir = catfile( $updir, $updir, "forest" );
} else {
  $forest_dir = catfile( $ENV{"CESR_SRC"}, "forest" );
}

if (-r catfile( $curdir, "tao", "code", "tao_struct.f90" )) {
  $tao_dir = catfile( $curdir, "tao" );
} elsif (-r catfile( $updir, "tao", "code", "tao_struct.f90")) {
  $tao_dir = catfile( $updir, "tao" );
} elsif (-r catfile( $updir, $updir, "tao", "code", "tao_struct.f90")) {
  $tao_dir = catfile( $updir, $updir, "tao" );
} else {
  $tao_dir = catfile( $ENV{"CESR_SRC"}, "tao" );
}

$str = @ARGV[0];
$str =~ s/\*/\\w\*/g;  # replace "*" by "\w*"

find(\&searchit, $bmad_dir);
find(\&searchit, $dcslib_dir);
find(\&searchit, $cesr_utils_dir);
find(\&searchit, $tao_dir);
## find(\&searchit, $recipes_dir);
## find(\&searchit, $forest_dir);

if ($found_one == 0) {print "Cannot match String! $str";}
print "\n";

#---------------------------------------------------------

sub searchit {

  if (/\#/) {return;}
  $file = $_;

# If a f90 file look in the file for a match.

  if ($file =~ /\.f90$/i) {

    $found_in_file = 0;

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

      # if a routine then look for match

      if (&routine_here) {

        # If a match then print info

        if ($routine_name =~ /^\s*$str[ \(\n]/i) {
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
        if ($1 =~ /^$str$/i) {
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

  if (/^$str\.f90$/ && $found_in_file == 0) {
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
