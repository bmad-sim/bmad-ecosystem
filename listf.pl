#!/usr/bin/perl

use File::Find;

$found_one = 0;

if (-d "./bmad")
  {$bmad_dir="./bmad";}
elsif (-d "../bmad")
  {$bmad_dir="../bmad";}
elsif (-d "../../bmad")
  {$bmad_dir="../../bmad";}
else
  {$bmad_dir=$ENV{"CESR_CVSSRC"}."/bmad";}

if (-d "./cesr_utils")
  {$cesr_utils_dir="./cesr_utils";}
elsif (-d "../cesr_utils")
  {$cesr_utils_dir="../cesr_utils";}
elsif (-d "../../cesr_utils")
  {$cesr_utils_dir="../../cesr_utils";}
else
  {$cesr_utils_dir=$ENV{"CESR_CVSSRC"}."/cesr_utils";}

if (-d "./dcslib")
  {$dcslib_dir="./dcslib";}
elsif (-d "../dcslib")
  {$dcslib_dir="../dcslib";}
elsif (-d "../../dcslib")
  {$dcslib_dir="../../dcslib";}
else
  {$dcslib_dir=$ENV{"CESR_CVSSRC"}."/dcslib";}
 
find(\&searchit, $bmad_dir);
find(\&searchit, $dcslib_dir);
find(\&searchit, $cesr_utils_dir);

if ($found_one == 0) {print "Cannot match String!";}
print "\n";

#---------------------------------------------------------

sub searchit {

  $file = $_;
  $str = @ARGV[0];
  $str =~ s/\*/\.\*/g;  # replace "*" by ".*"

# See if the file name matches.

  if (/^$str\.f90$/) {
    $found_one = 1;
    print "$File::Find::name\n";
  }

# If in the modules directory then look in the file for a match.

  $found_in_file = 0;

  if ($File::Find::name =~ m@modules/[^/]*\.f90$@ || $file =~ /mod\.f90$/) {

    open (F_IN, $file) || die ("Cannot open File: $_");
    while (<F_IN>) {

      if (/^ *interface$/i || /^ *interface /i) {   # skip interface blocks
        while (<F_IN>) {
          if (/^ * end interface/i) {last;}
        }
      }

      if (/^ *subroutine /i || /^ *recursive subroutine /i || 
            /^ *function /i || /^ *type /i || /^ *elemental subroutine /i ||
            /^ *real\(rp\) *function /i) {
        $name = $';              # strip off "subroutine
        $name =~ s/\(.*//;       # strip off "(..."
        if ($name =~ /^\s*$str\s*$/) {
          if ($found_in_file == 0) {print "Module: $File::Find::name\n";}
          $found_one = 1;
          $found_in_file = 1;
          print "     $_";
          if (/&\s*$/) {
            $_ = <F_IN>;
            print "   $_";
          }
        }
      }
    }
    close (F_IN);
  }
  
  $_ = $file;

}
