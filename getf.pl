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
 
$str = $field[$i];
$str =~ s/\*/\\w\*/g;   # replace "*" by "\w*" (any word characters)

#

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

# See if the file name matches.

  if ($file =~ /^$str\.f90$/i) {

    print "\n$File::Find::name\n";
    $found_one = 1;

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    $found = 0;

    while (<F_IN>) {
      if (/^!/) {
        print "$_";
        $found = 1; }
      elsif ($found == 1) {
        last;
      }
    }
  } elsif ($file =~ /\.f90$/) {
    $recording = 0;
    @comments = ();     

    $stat = open (F_IN, $file);
    if (! $stat) {
      print STDOUT "Cannot open File: $file\n"; 
      $_ = $file;
      return;
    }

    while (<F_IN>) {

      if (/^ *interface *$/i) {    # skip interface blocks
        while (<F_IN>) {
          if (/^ *end +interface/i) {last;}
        }
      }

      if (/^ *type *\(/i) {next;}   # skip "type (" constructs

      if (/^!\+/) {      # if match to "!+" start recording comment lines
        @comments = ($_);
        $recording = 1;

      } elsif ($recording == 1 && /^!\-/) {
        @comments = (@comments, $_);
        $recording = 0;

      } elsif ($recording == 1) {
        @comments = (@comments, $_)

      # match to type statement

      } elsif (/^ *type +$str[ \n]/i) {
        $found_one = 1;
        print "\n$File::Find::name\n";
        print $_;
        while (<F_IN>) {
          print $_;
          if (/^ *end type/i) {last;}
        }

      # match to subroutine, function, etc.

      } elsif (&routine_here) {
        $_ = $routine_name;     # strip off "subroutine"
        if (/^\s*$str[ |\(\n]/i) {
          $found_one = 1;
          print "\n$File::Find::name\n";
          foreach $line (@comments) {print $line;}
        }
        @comments = ();
        $recording = 0;

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
        if ($1 =~ /^$str$/i) {
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
