#!/usr/bin/perl -w

#############################################################
# it creates "data base file" which contain info where functions or
# subroutine are. Format of each line in the file:
# file.c : /dir/where/it/lives : name_of_function 
# ussage:  ./create_fdb.pl
# output: cesr_file.fdb
#############################################################
@tps_list=qw(h inc c cc cpp cxx f F f90 F90); #file extensions in use
open FILE1, ">".$ENV{"CESR_CONFIG"}."/cesr_file.fdb" || die "can not open file cesr_file.fdb";
#take all of the directories except CVS and make loop on them
foreach $dir (`ls -R \$CESR_SRC \$CESR_PKG | grep \':\'|grep -v .svn`) {
    chomp $dir;
    chop $dir;
    $ff=3;
    $ls_arg="";
    # create command line arg to see all files in the directory
    foreach $x (@tps_list){
        $ls_arg=$ls_arg."$dir/*\.$x ";
    }
    
    # Loop through all files in directory, throw out 'ls' complaints
    # about certain filetypes not being present (prevents tons of 
    # extraneous output.
    foreach $x  (`ls $ls_arg 2> /dev/null`){
        chomp $x;
        if ($x=~/\/(\w+\-*\w*\.\w+)/){$x=$1};#notdir filename
        if ($x=~/(?:\w+\-*\w*)\.(?:h|c|cc|cpp|cxx)/i){$ff=deal_c();}
        if ($x=~/(?:\w+\-*\w*)\.(?:inc|f|f90)/i){$ff=deal_f();}
        if ($x=~/(?:\w+\-*\w*)\.(?:h|c|cc|cpp|cxx|inc|f|f90)/i){
            if ($ff==0) {print FILE1 "$x : $dir : \n";}
        }
    }
} 
close FILE1 || die "couldn't close file: cesr_file.fdb";

#------------------------------------------------------------
# sub deal_c takes function definition from c-file and prints function name
sub deal_c{
    my ($file,$found_function);
    $found_function=0;
    $file="$dir/$x";
    open (FILE2,"<$file") || die "can not open file $file";
    $in_comment=0; # not in a comments field
  LINE:while (<FILE2>){
      $_ =~ s/\/\*.*\*\///g; #ignore inline comments of type /*...*/
      if ( $_=~ /\*\//) {
          $in_comment=0;
          $_=~s/.*\*\///g;
      }
      next LINE if $in_comment == 1;
      #search for function statement and printing function name
#      if ( $_ =~ /^(?:\s*\t*)(?:\w+\s+\t*){1,3}(\w+)\W/){
      if ( $_ =~ /^(?:\s*\t*)(?:(?:void|int|char|float|double|bool|wchar_t|signed|unsigned|long|short|extern|\"C\")\*?\s+\t*){1,3}\*?(\w+)(?:\[.*\])?(?:\s*\t*\(.*)/i){
              print FILE1 "$x : $dir : $1\n"; #basic definition
              $found_function=1;
      }
      if ( $_ =~ /^(?:\s*\t*)(?:(?:void|int|char|float|double|bool|wchar_t|signed|unsigned|long|short|extern|\"C\")\*?\s+\t*){1,3}\(\*?(\w+)(?:\[.*\])?\)(?:\s*\t*\(.*)/i ){
              print FILE1 "$x : $dir : $1\n";#definition type int* (*func_name)(...)
              $found_function=1;
      }
      if ( $_=~ /\/\*/) {
          $in_comment=1; #in a comments field
      }
  }
    close (FILE2) || die "can not close file $file";
    return $found_function;
}

# sub deal_f takes subroutine and function f90 definition from file and prints their name
sub deal_f{
    my ($file,$found_function);
    $found_function=0;
    $file="$dir/$x";
    open (FILE2,"<$file") || die "can not open file $file";
  LINE:while (<FILE2>){
     
      #search for function and subroutine statement and print their name
      if ( $_ =~ /^(?:\s*\t*)(?:(?:subroutine|function)\s+\t*)(\w+)\W/i){
              print FILE1 "$x : $dir : $1\n";
              $found_function=1;
      }
  }
    close (FILE2) || die "can not close file $file";
    return $found_function;
}
