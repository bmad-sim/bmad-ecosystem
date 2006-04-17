@echo off
rem the following two lines must be changed:
set LIBREL=bmad_dist
set BMAD_DRIVE=W:

set BMAD_REL=%BMAD_DRIVE%\%LIBREL%
set BMAD_MOD=%BMAD_REL%\modules
set BMAD_GMAKE=%BMAD_REL%\Gmake
SET BMAD_CONFIG=%BMAD_REL%\config
SET BMAD_LAYOUT=%BMAD_CONFIG%\layout
set BMAD_UTIL=%BMAD_REL%\util
set BMAD_inc=%BMAD_REL%\include
set BMAD_LAT=%BMAD_CONFIG%\lat
set BMAD_EXE=%BMAD_REL%\bin
set BMAD_LIB=%BMAD_REL%\lib
set PGPLOT_FONTS=%BMAD_REL%\PGPLOT

