$!########################################################
$!
$! Set CESR library logical variables for the new CESR
$! library scheme.  This is the VMS version which is 
$! based on the UNIX cesrenv.
$!
$! Mark Palmer  Oct. 24, 2001
$!
$!########################################################
$!
$!
$!----- Set up JOB versus global definition mode -----
$!
$ EXPORT = "define"
$! if p1.nes.""
$! then
$!   if (f$edit(p1,"UPCASE, TRIM").eqs."JOB")
$!   then
$!     write sys$output "In USER MODE"
$!     EXPORT = "define/job"
$!   endif
$! endif
$!
$!
$!-----  Outside packages  -----
$!
$ PBASE  = f$extract(0,f$locate("]",PACKAGES),PACKAGES)
$ EXPORT CESR_PKG "''PACKAGES'"
$ EXPORT PKG_LIB "''PBASE'.lib]"
$ EXPORT PKG_MOD "''PBASE'.modules]"
$!
$!
$!------  CESR libraries  ------
$!
$ CBASE  = f$extract(0,f$locate("]",CESR_),CESR_)
$ EXPORT CESR_CURRENT "''CBASE'.current]"
$ EXPORT CESR_DEVEL   "''CBASE'.devel]"
$!
$!
$!-- CESR Unix documentation --
$!
$ EXPORT CESR_DOC     "''CBASE'.doc]"
$!
$!
$!---  CESR CVS Repository DOESN'T EXIST on VMS  ---
$!
$! EXPORT CESR_CVSROOT "''CBASE'.cvsroot]"
$!
$!
$!----  CESR_LIB subtrees  ----
$!
$ EXPORT CESR_LIB      "''CBASE'.''CESRLIB']"
$ CLBASE  = "''CBASE'.''CESRLIB'"
$!
$ EXPORT CESR_CONFIG    "''CLBASE'.config]"
$ EXPORT CESR_CVSSRC    "''CLBASE'.cvssrc]"
$ EXPORT CESR_EXE       "''CLBASE'.bin]"
$ EXPORT CESR_GMAKE     "''CLBASE'.Gmake]"
$ EXPORT CESR_INC       "''CLBASE'.include]"
$ EXPORT CESR_MOD       "''CLBASE'.modules]"
$ EXPORT CESR_RUN       "''CLBASE'.run]"
$ EXPORT CESR_UTIL      "''CLBASE'.util]"
$!
$!
$!------  CESR constants  ------
$!
EXPORT CESR_CONST       "''CONST'"
$!
$!
$!-- BMAD configuration files --
$!
$ CONFBASE = "''CLBASE'.config"
$ EXPORT BMAD_LAT         "''CONFBASE'.bmad.lat]"
$ EXPORT BMAD_LAYOUT      "''CONFBASE'.bmad.layout]"
$!
$!
$! $Id$
$!
$! $Log$
$! Revision 1.2  2002/01/11 21:05:34  cesrulib
$! Bug fix in bmad_ logicals
$!
$! Revision 1.1  2002/01/03 22:21:44  cesrulib
$! Add scripts to set up CESR logicals on VMS.
$!
$! 
