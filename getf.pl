#!/usr/bin/perl

use File::Find;

$found_one = 0;

if (-d "./bmad/modules")
  {$bmad_dir="./bmad";}
elsif (-d "../bmad/modules")
  {$bmad_dir="../bmad";}
elsif (-d "../../bmad/modules")
  {$bmad_dir="../../bmad";}
else
  {$bmad_dir=$ENV{"CESR_CVSSRC"}."/bmad";}

if (-d "./cesr_utils/modules")
  {$cesr_utils_dir="./cesr_utils";}
elsif (-d "../cesr_utils/modules")
  {$cesr_utils_dir="../cesr_utils";}
elsif (-d "../../cesr_utils/modules")
  {$cesr_utils_dir="../../cesr_utils";}
else
  {$cesr_utils_dir=$ENV{"CESR_CVSSRC"}."/cesr_utils";}

if (-d "./dcslib/modules")
  {$dcslib_dir="./dcslib";}
elsif (-d "../dcslib/modules")
  {$dcslib_dir="../dcslib";}
elsif (-d "../../dcslib/modules")
  {$dcslib_dir="../../dcslib";}
else
  {$dcslib_dir=$ENV{"CESR_CVSSRC"}."/dcslib";}

if (-d "./recipes_f-90_LEPP")
  {$recipes_dir="./recipes_f-90_LEPP";}
elsif (-d "../recipes_f-90_LEPP")
  {$recipes_dir="../recipes_f-90_LEPP";}
elsif (-d "../../recipes_f-90_LEPP")
  {$recipes_dir="../../recipes_f-90_LEPP";}
else
  {$recipes_dir=$ENV{"CESR_CVSSRC"}."/recipes_f-90_LEPP";}

if (-d "./forest/basic")
  {$forest_dir="./forest/basic";}
elsif (-d "../forest/basic")
  {$forest_dir="../forest/basic";}
elsif (-d "../../forest/basic")
  {$forest_dir="../../forest/basic";}
else
  {$forest_dir=$ENV{"CESR_PKG"}."/forest/basic";}

if (-d "./tao/code")
  {$tao_dir="./tao";}
elsif (-d "../tao/code")
  {$tao_dir="../tao";}
elsif (-d "../../tao/code")
  {$tao_dir="../../tao";}
else
  {$tao_dir=$ENV{"CESR_CVSSRC"}."/tao";}

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

  $file = $_;

# See if the file name matches.

  if ($file =~ /^$str\.f90$/i) {

    print "\n$File::Find::name\n";
    $found_one = 1;
    open (F_IN, $_) || die ("Cannot open File: $_");

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
    open (F_IN, $file) || die ("Cannot open File: $_");
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
    open (F_IN, $file) || die ("Cannot open File: $_");

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

