#!/usr/bin/perl -w

#############################################################
# ussage:  ./find_file.pl file
#          ./find_file.pl file -I<dir1> -I<dir2> ... 
# it gives you location of the src file.
#############################################################

### taking argument from command line ###
if ($#ARGV == -1) {die "get_file.pl: no arguments";}
$argument ="$ARGV[0]";
$argument =~ s/\./\\\./;
$i=0;
while (<@ARGV>){
    if ($_=~ /-I(.*)/){
        $extra_dir_list[$i]=$1;
        $i=$i+1;
    }
}

### searching for the file in cesr_file.fdb
print "search $argument:\n";
open FILE1, "<".$ENV{"CESR_CONFIG"}."/cesr_file.fdb" || die "can not open file cesr_file.fdb ";

$i=0;
while (<FILE1>){
    if ($_=~ /^($argument\w*\.*\w*)(?:\s+\:\s+)(.+)(?:\s+\:\s+)(\w*)/){
        $file_list[$i]=$1;
        $release_dir_list[$i]=$2;
#       $func_list[$i]=$3;
        $local_dir_list[$i]=$release_dir_list[$i];
        $local_dir_list[$i]=~ s/.*\/cvssrc\//\.\.\//g;
        $local_dir_list[$i]=~ s/.*\/packages\//\.\.\//g;
        $i=$i+1;
    }
}
if ($#file_list == -1){  #exit if no file found in cesr_file.fdb 
    print "file not found\n";
    exit 0;
}
if (funic()==0){  #if found file is unic check local directories
    $i=0;         #and print them  
    if (-f "./$file_list[$i]"){
        print "Local File  : ./$file_list[$i]\n";
    }
    if (-f "$local_dir_list[$i]/$file_list[$i]"){
        print "Local File  : $local_dir_list[$i]/$file_list[$i]\n";
    }
    while (<@extra_dir_list>){
        if (-f "$_/$file_list[$i]") {
            print "Extra Local File  : $_/$file_list[$i]\n";
        }
    }
    print "Release File: $release_dir_list[$i]/$file_list[$i] \n";
    exit 0;
}

print "Multiple files found!:\n";
@sorted_list = sort @file_list;
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

while (<@shrink_list>) {#print all names of multiple files
    print "$_\n";
}

close FILE1 || die "couldn't close file: cesr_file.fdb";

#the following function checks if file_list contains only one unic name 
sub funic{ 
    my ($i,$key,$st);
    $key=0;
    $st=$file_list[0];
    while ( <@file_list> ) {
        if ( $st ne $_ ){ $key=1;}
    }
    return $key
    }
        
