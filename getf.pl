#!/usr/bin/perl

use File::Find;

$found_one = 0;

if(-d "./bmad")
{$bmad_dir="./bmad";}
elsif(-d "../bmad")
{$bmad_dir="../bmad";}
else
{$bmad_dir=$ENV{"CESR_CVSSRC"}."/bmad";}
if(-d "./cesr_utils")
{$cesr_utils_dir="./cesr_utils";}
elsif(-d "../cesr_utils")
{$cesr_utils_dir="../cesr_utils";}
else
{$cesr_utils_dir=$ENV{"CESR_CVSSRC"}."/cesr_utils";}
if(-d "./dcslib")
{$dcslib_dir="./dcslib";}
elsif(-d "../dcslib")
{$dcslib_dir="../dcslib";}
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
  $str =~ s/\*/\.\*/g;   # replace "*" by ".*"

# See if the file name matches.

  if (/^$str\.f90$/) {

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
  }

# If in the modules directory then look in the file for a match.

  if ($File::Find::name =~ m@modules/[^/]*\.f90$@) {
    $recording = 0;
    @comments = ();     
    open (F_IN, $file) || die ("Cannot open File: $_");
    while (<F_IN>) {

      if (/^ *interface$/i || /^ *interface /i) {    # skip interface blocks
        while (<F_IN>) {
          if (/^ * end interface/i) {last;}
        }
      }

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
             /^ *real\(rp\) *function /i) {
        $_ = $';     # strip off "subroutine"
        s/\(.*//;    # strip off "(..."
        if (/^\s*$str\s*$/) {
          $found_one = 1;
          print "\n$File::Find::name\n";
          foreach $line (@comments) {print $line;}
        }
        @comments = ();
        $recording = 0;
      }
      elsif (/^ *type *$str/) {
        $found_one = 1;
        print "\n$File::Find::name\n";
        print $_;
        while (<F_IN>) {
          print $_;
          if (/^ * end type/) {last;}
        }
      }
    }
    close (F_IN);
  }
  

  $_ = $file;

}
