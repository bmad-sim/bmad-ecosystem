#!/usr/bin/perl -w

###########################################################
# ussage: ./f90_mod_rules.perl file0 file1 file2 ...,
# Search file for modules and print modules names 
###########################################################

###########################################################
######### taking arguments from comand line #######
if ($#ARGV == -1) {die "f90_mod_list.perl: no arguments";}

######### searching file for modules and printing names ###
while ($file=<@ARGV>) {
open (FILE1,"$file")|| die "can not open file $file";
$x="";
while (<FILE1>){
    if ( $_ =~ /^(?:\s*\t*)module(?:\s+\t*)(\w+).*/i ){
        $x=$1;
        if ($x !~ /procedure/i){
        print "$x ";
    }
    }
}
close(FILE1) || die "couldn't close $file";
}














