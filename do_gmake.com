$!################################################################
$!
$! Run gmake in a specified sub-directory
$! This file allows us to get around the issue that shell, 
$! wildcard, and vpath functionality does not (yet?) work 
$! properly in the VMS version of GNU Make.
$! 
$! Usage:  @[cesrulib.com]vms_sources <mode> <dir> <variable_def>
$!
$!         <source_dest> =   OPN  start a new file
$!                           VAR  append a variable definition
$!                           OLB  add sources being compiled into
$!                                an object library
$!                           OBJ  add sources being compiled into
$!                                individual object files
$!                           MOD  add module sources
$!                           CFG  add configuration files
$!
$! Author:  Mark Palmer  Jan 3, 2002
$!
$!################################################################
$!
$ set noverify
$!
$!-----------------------
$! Check for valid input
$!-----------------------
$!
$ tp1 = 0
$ tp2 = 0
$ tp3 = 0
$ if ("''p1'" .eqs. "") then tp1 = 1
$ if ("''p2'" .eqs. "") then tp2 = 1
$ if ("''p3'" .eqs. "") then tp3 = 1
$ if (tp1.or.tp2)
$ then 
$   write sys$output "Incorrect inputs specified..."
$   exit
$ endif
$!
$!-------------------------------------
$! Move to the specified directory
$!-------------------------------------
$!
$ set def 'p1'
$!
$!-------------------------------------
$! Run gmake
$!-------------------------------------
$!
$ if (tp3) 
$ then
$   gmake 'p3' -f 'p2'
$ else
$   gmake -f 'p2'
$ endif
$!
$ exit
$!
$!################################################################
$!
$! $Id$
$!
$! $Log$
$! Revision 1.1  2002/01/09 16:29:57  cesrulib
$! First pass at VMS library build infrastructure
$!
$!
$!################################################################


