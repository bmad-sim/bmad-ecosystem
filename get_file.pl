#!/usr/bin/perl -w

#############################################################
# ussage:  ./get_file.pl  [function_name]
#          ./get_file.pl  [subs_name]
# output:  dir/filename
# it gives you info where function or subroutine is defined or
# declared. 
#############################################################

### taking argument from command line ###
if ($#ARGV == -1) {die "get_file.pl: no arguments";}
if ($#ARGV == 1) {die "get_file.pl: more then one argument";}
$argument ="$ARGV[0]";

### searching for the function
print "search for $argument:\n";
open FILE1, "<".$ENV{"CESR_CONFIG"}."/cesr_file.fdb" || die "can not open file cesr_file.fdb ";
while (<FILE1>){
    if ($_=~ /^(\w+\.\w+)(?:\s+\:\s+)(.+)(?:\s+\:\s+)$argument/){
        print "$2/$1\n";
    }
}
close FILE1 || die "couldn't close file: cesr_file.fdb";


