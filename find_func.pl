#!/usr/bin/perl -w

#############################################################
# ussage:  ./find_func.pl function
#          ./find_func.pl function -I<dir1> -I<dir2> ... 
# it gives location of the function.
#############################################################

### taking argument from command line ###
if ($#ARGV == -1) {die "get_func.pl: no arguments";}
$argument ="$ARGV[0]";
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
    if ($_=~ /^(\w+\.\w+)(?:\s+\:\s+)(.+)(?:\s+\:\s+)($argument\w*)/){
        $file_list[$i]=$1;
        $release_dir_list[$i]=$2;
        $func_list[$i]=$3;
        $local_dir_list[$i]=$release_dir_list[$i];
        $local_dir_list[$i]=~ s/.*\/cvssrc\//\.\.\//g;
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
    while (<@local_dir_list>){
        if (-f "$_/$file_list[$i]"){
            print "Local File  : $_/$file_list[$i]\n";
            $i=$i+1;
        }
    }
    $i=0;
    while (<@extra_dir_list>){
        if (-f "$_/$file_list[$i]") {
            print "Extra Local File  : $_/$file_list[$i]\n";
            $i=$i+1;
        }
    }
    $i=0;
    while (<@release_dir_list>){
        print "Release File: $_/$file_list[$i] \n";
        $i=$i+1;
    }
    exit 0;
}

print "Multiple function names found!:\n";
@sorted_list = sort @func_list;
$i=0;
$sl=$sorted_list[0]; #shrink sorted list (remove repeated file_names).
$shrink_list[0]=$sl;
while (<@sorted_list>){
    if ($_ ne $sl){
        $sl=$_;
        $i=$i+1;
        $shrink_list[$i]=$sl;
    }

}

while (<@shrink_list>) {#print all multiple function names
    print "$_\n";
}
print "--------------\n";
while (<@func_list>) {#print all multiple function names
    print "$_\n";
}

close FILE1 || die "couldn't close file: cesr_file.fdb";

#the following function checks if func_list contains only one unic name 
sub funic{ 
    my ($i,$key,$st);
    $key=0;
    $st=$func_list[0];
    while ( <@func_list> ) {
        if ( $st ne $_ ){ $key=1;}
    }
    return $key
    }
        
