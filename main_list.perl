#!/usr/bin/perl -w

###########################################################
# ussage: ./main_list.perl file0 file1 file2 ...,
# Search files for main programs and print file_base_names.
# If there is multiple definitions of the "main" in a file it
# will print filename only once 
###########################################################

###########################################################
######### checking arguments from comand line #######
if ($#ARGV == -1) {die "f90_mod_list.perl: no arguments";}

## searching file for main (program) and 
##         printing file_base_names
##
while ($file=<@ARGV>) {
open (FILE1,"$file")|| die "can not open file $file";
$in_comment=0; # not in a comments field
LINE:while (<FILE1>){
    $_ =~ s/\/\*.*\*\///g; #ignore inline comments of type /*...*/
    if ( $_=~ /\*\//) {
        $in_comment=0;
        $_=~s/.*\*\///g;
    }
    next LINE if $in_comment == 1;
    #search for "program" statement and printing file_base_name
    if ($_ =~ /^(?:\s*\t*)program/i){
        if ($file =~ /(\w+)\.f(\w*)\Z/i){ #file_base_name
        print "$1 ";
        last LINE;
    }
    }
    #search for "main" statement and printing file_base_name
    if ( $_ =~ /^(?:\s*\t*)(?:\w*\s+\t*){0,2}main\W/i){
        if ($file =~ /(\w+)\.(c|cc|cpp|cxx)\Z/i){ #file_base_name
        print "$1 ";
        last LINE;
    }
    }
    if ( $_=~ /\/\*/) {
        $in_comment=1; #in a comments field
    }
}
close(FILE1) || die "couldn't close $file";
}














