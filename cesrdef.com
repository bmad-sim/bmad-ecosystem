$!---- cesrdef.com ----
$!
$! Define cesr specific logicals;
$! to be used on VMS.
$!
$! Shamelessly stolen from the CLEO cleodef script
$! M. Palmer
$!
$! Usage in .login:
$!   source /home/cesrulib/bin/cesrdef
$! or switch to other versions of libraries and constants:
$!   setenv CESRLIB devel; source /usr/local/bin/cesrdef
$! or define aliases in .cshrc, e.g.:
$!   alias cesrdef "source /home/cesrulib/bin/cesrdef"
$!   alias current "setenv CESRLIB current ; cesrdef"
$!   alias devel   "setenv CESRLIB devel ; cesrdef"
$!--------------------------------------------------------
$!
$!--------------------------------------------------------
$!
$!       Root names
$!
$!--------------------------------------------------------
$!
$! Check whether CONST is defined
$!
$ set message/noseverity/nofacility/noident/notext
$ sho sym CONST
$ defined = $status
$ set message/severity/facility/ident/text
$!
$ if .not.defined
$ then
$  CONST = ""
$ endif
$!
$! Now set CONST
$!
$ if (CONST.nes."") 
$ then
$!  write sys$output "CONST    =  "+CONST
$ else
$  CONST = "U:[cesr.constants]"
$  write sys$output "  CONST =  ""''CONST'"" "
$ endif
$!
$!
$! Check whether CESR_ is defined
$!
$ set message/noseverity/nofacility/noident/notext
$ sho sym CESR_
$ defined = $status
$ set message/severity/facility/ident/text
$!
$ if .not.defined
$ then
$  CESR_ = ""
$ endif
$!
$! Now set CESR_
$!
$ if (CESR_.nes."") 
$ then
$!  write sys$output "CESR_    =  "+CESR_
$ else
$  CESR_ = "U:[cesrulib.cesr_libs]"
$  write sys$output "  CESR_ =  ""''CESR_'"" "
$ endif
$!
$!
$! Check whether PACKAGES is defined
$!
$ set message/noseverity/nofacility/noident/notext
$ sho sym PACKAGES
$ defined = $status
$ set message/severity/facility/ident/text
$!
$ if .not.defined
$ then
$  PACKAGES = ""
$ endif
$!
$! Now set PACKAGES
$!
$ if (PACKAGES.nes."") 
$ then
$!  write sys$output "PACKAGES =  "+PACKAGES
$ else
$  PACKAGES = "U:[cesrulib.cesr_libs.packages]"
$  write sys$output "  PACKAGES =  ""''PACKAGES'"" "
$ endif
$!
$!
$!       Switchable version
$!
$!
$! Check whether CESRLIB is defined
$!
$ set message/noseverity/nofacility/noident/notext
$ sho sym CESRLIB
$ defined = $status
$ set message/severity/facility/ident/text
$!
$ if .not.defined
$ then
$  CESRLIB = ""
$ endif

$!
$! Now set CESRLIB
$!
$ if (CESRLIB.nes."") 
$ then
$  write sys$output "CESRLIB  =  "+CESRLIB
$ else
$  CESRLIB = "current"
$  write sys$output "CESRLIB  =  "+CESRLIB
$ endif
$!
$!
$ CESRENV = "U:[CESRULIB.BIN]CESRENV.COM"
$ write sys$output CESRENV
$ if (f$search(CESRENV).nes."")
$ then
$  @'CESRENV'
$ else
$  write sys$output CESRENV+" not found."
$ endif
$!
$! $Id$
$!
$! $Log$
$! Revision 1.1  2002/01/03 22:21:44  cesrulib
$! Add scripts to set up CESR logicals on VMS.
$!
$!




