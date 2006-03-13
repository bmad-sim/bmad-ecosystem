#!/usr/bin/perl -w


###########################################################
# ussage:  ./f90_mod_rules.perl file0 file1 file2 ... flagg,
# example: ./f90_mod_rules.perl file0 file1 file2 ... u,
#          ./f90_mod_rules.perl file0 file1 file2 ... l,
# flagg = u  (upper), or flagg = anything (lower); 
# Search file for modules and print modules names in upper
# or lower case. 
# 
###########################################################

###########################################################
######### taking arguments from comand line #######
if ($#ARGV == -1) {die "f90_mod_list.perl: no arguments";}
if ($#ARGV == 0) {die "f90_mod_list.perl: just one argument";}
$flagg = pop(@ARGV);
$case ="lower";
if ($flagg eq "u"){$case="upper"};

######### searching file for modules and printing names ###
while ($file=<@ARGV>) {
open (FILE1,"$file")|| die "can not open file $file";
$x="";
while (<FILE1>){
    if ( $_ =~ /^(?:\s*\t*)module(?:\s+\t*)(\w+).*/i ){
        $x=$1;
        if ($x !~ /procedure/i){
           if ($case eq "upper"){print "\U$x\E "} 
           else{ print "\L$x\E "}
    }
    }
}
close(FILE1) || die "couldn't close $file";
}














