#!/usr/bin/perl -w

###########################################################
# ussage: ./f90_mod_rules.perl RULE OUT_DIR FILE,
# Search file for modules and gives output as a gmake rule 
# to compile modules from the source file.
#    
###########################################################

###########################################################
######### taking arguments from comand line #######

if ($#ARGV != 2) {die "f90_mod_rules.perl: wrong arguments";}
$RULE=$ARGV[0];
$OUT_DIR=$ARGV[1];
$absfilename=$ARGV[2];
$filename = $absfilename;
if ($absfilename =~ /(\w+\.{0,1}\w+)\Z/){$filename="$1"}
######### searching file for modules  ###############
open (FILE1,"$absfilename")|| die "can not open file $absfilename";
$x="";
while (<FILE1>){
    if ( $_ =~ /^(?:\s*\t*)module(?:\s+\t*)(\w+).*/i ){
        $x=$1;
        if ($x !~ /procedure/i){
        print "$OUT_DIR/\L$x\E.mod: $filename\n";
        print "\t$RULE\n";
    }
    }
}
close(FILE1) || die "couldn't close $absfilename";















