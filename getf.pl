#!/usr/bin/perl
# modify to use File::Spec for transparent use on VMS
# 2005.09.28 mjf

use File::Find;
use File::Spec::Functions;

# Following two lines should be uncommented on VMS to check the master_list:
#my @args = ("type_header", @ARGV);
#system(@args) == 0 or warn "system @args failed: $!\n";

$found_one = 0;

# -d is a file test operator which checks if the file is a directory

my $curdir = curdir();
my $updir = updir();

if (-d catdir( $curdir, "bmad", "modules" )) {
  $bmad_dir = catdir( $curdir, "bmad" );
} elsif (-d catdir( $updir, "bmad", "modules")) {
  $bmad_dir = catdir( $updir, "bmad" );
} elsif (-d catdir( $updir, $updir, "bmad", "modules")) {
  $bmad_dir = catdir( $updir,$updir, "bmad" );
} else {
  $bmad_dir = catdir( $ENV{"CESR_CVSSRC"}, "bmad" );
}

if (-d catdir( $curdir, "cesr_utils", "modules" )) {
  $cesr_utils_dir = catdir( $curdir, "cesr_utils" );
} elsif (-d catdir( $updir, "cesr_utils", "modules")) {
  $cesr_utils_dir = catdir( $updir, "cesr_utils" );
} elsif (-d catdir( $updir, $updir, "cesr_utils", "modules")) {
  $cesr_utils_dir = catdir( $updir,$updir, "cesr_utils" );
} else {
  $cesr_utils_dir = catdir( $ENV{"CESR_CVSSRC"}, "cesr_utils" );
}


if (-d catdir( $curdir, "dcslib", "modules" )) {
  $dcslib_dir = catdir( $curdir, "dcslib" );
} elsif (-d catdir( $updir, "dcslib", "modules")) {
  $dcslib_dir = catdir( $updir, "dcslib" );
} elsif (-d catdir( $updir, $updir, "dcslib", "modules")) {
  $dcslib_dir = catdir( $updir,$updir, "dcslib" );
} else {
  $dcslib_dir = catdir( $ENV{"CESR_CVSSRC"}, "dcslib" );
}


if (-d catdir( $curdir, "recipes_f-90_LEPP" )) {
  $recipes_dir = catdir( $curdir, "recipes_f-90_LEPP" );
} elsif (-d catdir( $updir, "recipes_f-90_LEPP")) {
  $recipes_dir = catdir( $updir, "recipes_f-90_LEPP" );
} elsif (-d catdir( $updir, $updir, "recipes_f-90_LEPP")) {
  $recipes_dir = catdir( $updir, $updir, "recipes_f-90_LEPP" );
} else {
  $recipes_dir = catdir( $ENV{"CESR_CVSSRC"}, "recipes_f-90_LEPP" );
}


if (-d catdir( $curdir, "forest", "basic" )) {
  $forest_dir = catdir( $curdir, "forest", "basic" );
} elsif (-d catdir( $updir, "forest", "basic")) {
  $forest_dir = catdir( $updir, "forest", "basic" );
} elsif (-d catdir( $updir, $updir, "forest", "basic")) {
  $forest_dir = catdir( $updir,$updir, "forest", "basic" );
} else {
  $forest_dir = catdir( $ENV{"CESR_CVSSRC"}, "forest" );
}

if (-d catdir( $curdir, "tao", "code" )) {
  $tao_dir = catdir( $curdir, "tao" );
} elsif (-d catdir( $updir, "tao", "code")) {
  $tao_dir = catdir( $updir, "tao" );
} elsif (-d catdir( $updir, $updir, "tao", "code")) {
  $tao_dir = catdir( $updir,$updir, "tao" );
} else {
  $tao_dir = catdir( $ENV{"CESR_CVSSRC"}, "tao" );
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

# find(\&searchit, $recipes_dir);
# find(\&searchit, $forest_dir);

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
      }
      elsif ($recording == 1 && /^!\-/) {
        @comments = (@comments, $_);
        $recording = 0;
      }
      elsif ($recording == 1) {
        @comments = (@comments, $_)
      }
      elsif (/^ *subroutine /i || /^ *recursive subroutine /i || 
             /^ *function /i || /^ *elemental subroutine /i ||
             /^ *real\(rp\) *function /i || /^ *interface /i) {
        $_ = $';     # strip off "subroutine"
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

      # match to type statement

      elsif (/^ *type +$str[ \n]/i) {
        $found_one = 1;
        print "\n$File::Find::name\n";
        print $_;
        while (<F_IN>) {
          print $_;
          if (/^ *end type/i) {last;}
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



