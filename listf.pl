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

$str = @ARGV[0];
$str =~ s/\*/\\w\*/g;  # replace "*" by "\w*"

find(\&searchit, $bmad_dir);
find(\&searchit, $dcslib_dir);
find(\&searchit, $cesr_utils_dir);
find(\&searchit, $tao_dir);

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

      if (/^ *subroutine /i || /^ *recursive subroutine /i || 
            /^ *function /i || /^ *type /i || /^ *elemental subroutine /i ||
            /^ *real\(rp\) *function /i || /^ *interface /i) {
        $name = $';              # strip off "routine" string

        # If a match then print info

        if ($name =~ /^\s*$str[ \(\n]/i) {
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
            $_ = $';  
            if (/^ *subroutine/i || /^ *function/i ||
                /^ *type/i || /^ *interface/i) {$count = $count - 1;}
          }
          elsif (/^ *subroutine /i || /^ *recursive subroutine /i || 
                /^ *function /i || /^ *elemental subroutine /i ||
                /^ *real\(rp\) *function /i || /^ *type/i) {
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
