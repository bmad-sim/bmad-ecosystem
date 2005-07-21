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
$argument =~ s/\*/\\w\*\\.\*\\w\*/;
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
    if ($_=~ /^($argument)(?:\s+\:\s+)(.+)(?:\s+\:\s+)(\w*)/){
        $file_list[$i]=$1;
        $release_dir_list[$i]=$2;
#       $func_list[$i]=$3;
        $local_dir_list[$i]=$release_dir_list[$i];
        $local_dir_list[$i]=~ s/.*\/cvssrc\//\.\.\//g;
        $local_dir_list[$i]=~ s/.*\/packages\//\.\.\//g;
        $release_file_list[$i]="$release_dir_list[$i]/$file_list[$i]";
        $local_file_list[$i]="$local_dir_list[$i]/$file_list[$i]";
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

    @tmp_list = @local_file_list;shrink();@shrink_list=@shrink_tmp_list;
    foreach $i (@shrink_list){
        $place_g="Local File";
        if (-f "$i"){
            print "$place_g: $i \n";
        }
    }
    foreach $j (<@extra_dir_list>){
        if (-f "$j/$file_list[$i]") {
            print "Extra Local File  : $j/$file_list[$i]\n";
        }
    }
    @tmp_list = @release_file_list;shrink();@shrink_list=@shrink_tmp_list;
    foreach $i (@shrink_list){
        $place_g="Release File";
        if ($i=~/\/packages\//){$place_g="Packages File"};
        print "$place_g: $i \n";
    }
    exit 0;
}

print "Multiple files found!:\n";

@tmp_list = @file_list;shrink();@shrink_list=@shrink_tmp_list;
foreach $j (@shrink_list) {#print all names of multiple files
    print "$j\n";
}

close FILE1 || die "couldn't close file: cesr_file.fdb";

#the following function checks if file_list contains only one unic name 
sub funic{ 
    my ($i,$key,$st);
    $key=0;
    $st=$file_list[0];
    foreach $i ( @file_list ) {
        if ( $st ne $i ){ $key=1;}
    }
    return $key
    }
#the following function shrink array so there are no repeated elements        
sub shrink{
    my ($i,$sl,@sorted_list,$j);
    @sorted_list = sort @tmp_list;
    $i=0;
    $sl=$sorted_list[0]; #shrink sorted list (remove repeated file_names).
    $shrink_tmp_list[0]=$sl;
    foreach $j (@sorted_list){
        if ($j ne $sl){
            $sl=$j;
            $i=$i+1;
            $shrink_tmp_list[$i]=$sl;
        }

    }
}
