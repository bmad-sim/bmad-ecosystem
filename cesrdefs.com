$!---- cesrdefs.com ----
$!
$! Define CESR specific logical variables;
$! to be used on VMS.
$!
$! Wrapper com file that executes cesrdef.com
$!
$!
@[cesrulib.bin]cesrdef
$!
$!
$! $Id$
$!
$! $Log$
$! Revision 1.1  2002/01/03 22:21:44  cesrulib
$! Add scripts to set up CESR logicals on VMS.
$!
$!
