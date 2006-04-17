@echo off
rem the following two lines must be changed:
set LIBREL=cesr_2006_0411_d
set CESRDRIVE=W:

set CESR_DOC=%CESRDRIVE%\doc
set CESR_CONST=%CESRDRIVE%\constants
set CESR_CVSROOT=%CESRDRIVE%\cvsroot
set CESR_BASE=%CESRDRIVE%\CYGWIN_NT_i686
set CESR_REL=%CESR_BASE%\%LIBREL%
set CESR_DEVEL=%CESR_BASE%\devel
set CESR_CURRENT=%CESR_BASE%\current
set CESR_CVSSRC=%CESR_REL%\cvssrc
set CESR_INC=%CESR_CVSSRC%\include
set CESR_RUN=%CESR_REL%\run
set CESR_UTIL=%CESR_REL%\util
set CESR_LIB=%CESR_REL%\lib
set CESR_GMAKE=%CESR_REL%\Gmake
set CESR_PKG=%CESR_REL%\packages
set CESR_MOD=%CESR_REL%\modules
set CESR_CONFIG=%CESR_REL%\config
set CESR_EXE=%CESR_REL%\bin
set PKG_MOD=%CESR_PKG%\modules
set PKG_LIB=%CESR_PKG%\lib
set PKG_INC=%CESR_PKG%\include
set TOAD_DIR=%CESR_CVSSRC%\toad


SET BMAD_CONFIG=%CESR_CONFIG%\bmad
SET BMAD_LAT=%BMAD_CONFIG%\lat
SET BMAD_LAYOUT=%BMAD_CONFIG%\layout
set TMPDIR=%TMP%
