#!/usr/bin/perl -w

#############################################################
# ussage:  ./get_function.pl file
# output:  "file":  dir/file --- "fun/sub": function/subroutine 
# it gives you info what function or subroutine is defined or
# declared in the file.
#############################################################

### taking argument from command line ###
if ($#ARGV == -1) {die "get_file.pl: no arguments";}
if ($#ARGV == 1) {die "get_file.pl: more then one argument";}
$argument ="$ARGV[0]";

### searching for the file
print "search $argument:\n";
open FILE1, "<".$ENV{"CESR_CONFIG"}."/cesr_file.fdb" || die "can not open file cesr_file.fdb ";
while (<FILE1>){
    if ($_=~ /^($argument\w*\.*\w*)(?:\s+\:\s+)(.+)(?:\s+\:\s+)(\w+)/){
        print "file:  $2/$1 --- fun/sub:  $3\n";
    }
}
close FILE1 || die "couldn't close file: cesr_file.fdb";
