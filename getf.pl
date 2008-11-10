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

if (-r catfile( $curdir, "pecklib", "code", "butout.f90" )) {
  $pecklib_dir = catfile( $curdir, "pecklib" );
} elsif (-r catfile( $updir, "pecklib", "code", "butout.f90")) {
  $pecklib_dir = catfile( $updir, "pecklib" );
} elsif (-r catfile( $updir, $updir, "pecklib", "code", "butout.f90")) {
  $pecklib_dir = catfile( $updir, $updir, "pecklib" );
} else {
  $pecklib_dir = catfile( $ENV{"CESR_SRC"}, "pecklib" );
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

# Look for arguments

$extra = 0;

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

  last;

}

#
 
$match_str = $field[$i];
$match_str =~ s/\*/\\w\*/g;   # replace "*" by "\w*" (any word characters)

#

find(\&searchit, $bmad_dir);
find(\&searchit, $dcslib_dir);
find(\&searchit, $cesr_utils_dir);
find(\&searchit, $tao_dir);
find(\&searchit, $pecklib_dir);
## find(\&searchit, $recipes_dir);
## find(\&searchit, $forest_dir);

if ($found_one == 0) {print "Cannot match String! $match_str";}
print "\n";

#---------------------------------------------------------

sub searchit {

  # Do not search this directory.
  if ($File::Find::name =~ /pecklib\/original_from_vms/) {return;}  

  if (/\#/) {return;}
  $file = $_;
  $matches_file_name = 0;

  # Check contents of .f90 file.

  if ($file =~ /\.f90$/) {
    $n_com = 0;
    @comments = ();     

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    while (<F_IN>) {

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
