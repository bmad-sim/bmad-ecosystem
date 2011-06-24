#!/usr/bin/perl
# modify to use File::Spec for transparent use on VMS
# 2006.01.18 mjf

use File::Find;
use File::Spec::Functions;

# this is so perl will find getf_common.pl in the same directory as listf.pl

use FindBin ();
use lib "$FindBin::Bin";

# Now load getf_common.pl and call setup_dirs routine

require "getf_common.pl";

&setup_dirs;

# Following two lines should be uncommented on VMS to check the master_list:

#my @args = ("list_list", @ARGV);
#system(@args) == 0 or warn "system @args failed: $!\n";

&search_it ("listf");

#---------------------------------------------------------

sub searchit {

  # Do not search this directory.
  if ($File::Find::name =~ /cesr_utils\/original_from_vms/) {return;}  

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

      if (&routine_here(1)) {

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
          elsif (&routine_here(1)) {
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

