#!/usr/bin/perl -w

#############################################################
# ussage:  ./find_func.pl function
#          ./find_func.pl function -I<dir1> -I<dir2> ... 
# it gives location of the function.
#############################################################

### taking argument from command line ###
if ($#ARGV == -1) {die "get_func.pl: no arguments";}
$argument ="$ARGV[0]";
$argument =~ s/\*/\\w\*/;
$i=0;
while (<@ARGV>){
    if ($_=~ /-I(.*)/){
        $extra_dir_list[$i]=$1;
        $i=$i+1;
    }
}

### searching for the file in cesr_file.fdb
open FILE1, "<".$ENV{"CESR_CONFIG"}."/cesr_file.fdb" || die "can not open file cesr_file.fdb ";

$i=0;
while (<FILE1>){
    if ($_=~ /^(\w+\.\w+)(?:\s+\:\s+)(.+)(?:\s+\:\s+)($argument)/){
        $file_list[$i]=$1;
        $release_dir_list[$i]=$2;
        $func_list[$i]=$3;
        $local_dir_list[$i]=$release_dir_list[$i];
        $local_dir_list[$i]=~ s/.*\/cvssrc\//\.\.\//g;
        $local_dir_list[$i]=~ s/.*\/packages\//\.\.\//g;
        $i=$i+1;
    }
}
if ($#func_list == -1){  #exit if no function name found in cesr_file.fdb 
    print "function not found\n";
    exit 0;
}
if (funic()==0){  #if found func_list contains unic name check local directories
    print "function $func_list[0]:\n"; 
    $i=0;         #and print them  
    if (-f "./$file_list[$i]"){
        print "Local File  : ./$file_list[$i]\n";
    }
    foreach $j (@local_dir_list){
        if (-f "$j/$file_list[$i]"){
            print "Local File  : $j/$file_list[$i]\n";
            $i=$i+1;
        }
    }
    $i=0;
    foreach $j (@extra_dir_list){
        if (-f "$j/$file_list[$i]") {
            print "Extra Local File  : $j/$file_list[$i]\n";
            $i=$i+1;
        }
    }
    $i=0;
    foreach $j (@release_dir_list){
        $place_g="Release File";
        if ($j=~/\/packages\//){$place_g="Packages File"};
        print "$place_g: $j/$file_list[$i] \n";
        $i=$i+1;
    }
    exit 0;
}

print "Multiple function names found!:\n";
@sorted_list = sort @func_list;
$i=0;
$sl=$sorted_list[0]; #shrink sorted list (remove repeated file_names).
$shrink_list[0]=$sl;
foreach $j (@sorted_list){
    if ($j ne $sl){
        $sl=$j;
        $i=$i+1;
        $shrink_list[$i]=$sl;
    }

}

foreach $j (@shrink_list) {#print all multiple function names
    print "$j\n";
}

close FILE1 || die "couldn't close file: cesr_file.fdb";

#the following function checks if func_list contains only one unic name 
sub funic{ 
    my ($i,$key,$st);
    $key=0;
    $st=$func_list[0];
    foreach $i ( @func_list ) {
        if ( $st ne $i ){ $key=1;}
    }
    return $key
    }
        
