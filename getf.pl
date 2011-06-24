#!/usr/bin/perl
# modify to use File::Spec for transparent use on VMS
# 2005.09.28 mjf

use File::Find;
use File::Spec::Functions;


# this is so perl will find getf_common.pl in the same directory as getf.pl

use FindBin ();
use lib "$FindBin::Bin";

# Now load getf_common.pl and call setup_dirs routine

require "getf_common.pl";

&setup_dirs;

# Following two lines should be uncommented on VMS to check the master_list:

#my @args = ("type_header", @ARGV);
#system(@args) == 0 or warn "system @args failed: $!\n";

$found_one = 0;

&search_it ("getf");

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

      } elsif (&routine_here(0)) {
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
          elsif (&routine_here(0)) {
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
      if (&routine_here(0)) {last;}
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
