$!################################################################
$!
$! Check for the existence of an object library file and create 
$! as necessary
$!
$! Mark Palmer  Jan 3, 2002
$!
$!################################################################
$ set noverify
$!
$! Check for valid input
$!
$ if ("''p1'" .eqs. "") 
$ then 
$   write sys$output "No library specified..."
$   exit
$ endif
$!
$! Check for library and create if not present.  Note that this 
$! com file presumes that higher directories in the hierarchy are
$! already present!
$!
$ if f$search(p1) .eqs. ""
$ then
$   write sys$output "Creating object library:  ", p1
$   lib/create 'p1'
$ endif
$!
$! $Id$
$!
$! $Log$
$! Revision 1.1  2002/01/09 16:29:57  cesrulib
$! First pass at VMS library build infrastructure
$!
$!
