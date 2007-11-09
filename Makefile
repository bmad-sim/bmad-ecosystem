#=============================================================================
#
# File:     bmad/Makefile
#
# Purpose:  Makefile to build the CESR BMAD Library
#
# Author:   M. Palmer   7/24/00
#
# Acknowledgements:  Simon Patton's CLEO Makefile infrastructure
#                    Dave Sagan's descrip.mms 
#                    Mike Marsh and his makefiles
#
# Modifications:  MAP - 7/24/00-9/23/01 - Many bug fixes and improvements
#
#=============================================================================
#
# Revision 1.11  2005/01/21 13:53:13  palmer
# Update bmad makefiles for addition of test_job area with sources for a
# standard test program.
#
# Revision 1.10  2004/08/26 15:03:05  dcs
# Added cpp_interface directory.
#
# Revision 1.9  2004/04/19 23:07:06  dcs
# Added bmad/cesr directory.
#
# Revision 1.8  2003/11/10 14:05:08  palmer
# Update for BMAD distributions.
#
# Revision 1.7  2003/04/21 17:24:14  cesrulib
# new makefile
#
# Revision 1.6  2002/06/13 14:53:29  dcs
# Interfaced with FPP/PTC
#
# Revision 1.5  2002/02/23 20:31:35  dcs
# Double/Single Real toggle added
#
# Revision 1.4  2001/10/23 21:34:14  palmer
# Update to work with M.tail infrastructure.
#
# Revision 1.3  2001/10/05 12:26:33  palmer
# More bug fixes.
#
# Revision 1.2  2001/10/05 04:08:45  palmer
# Bug fixes.
#
# Revision 1.1  2001/10/05 02:43:28  palmer
# Install the Makefile for the BMAD library.
#
#

#-----------------------------------------------------------------------------
# Specify default source directories for code to be compiled and put into
# libraries, for module files, for code to be compiled into explicit object
# files, and initialization files to be copied to the config directory tree.
# Also specify lists of LOCAL directories to search for include files and 
# compiled module files (NOTE:  the directories $(CESR_INC) and $(CESR_MOD)
# are automatically included in this list by M.tail).  
#
# Multiple directories should be separated by spaces.
#
#
# WARNING:  DO NOT put a trailing SLASH on directory names
#
# WARNING:  EMPTY lists should have NO EXTRA SPACES after the :=
#           It is better to simply delete them.  The ones shown here 
#           are intended to show the full range of options available.
#
#
# CODE_SRC_DIRS   - Put a list of code sub-directories here
# MOD_SRC_DIRS    - Put a list of F90 module source sub-directories here
# OBJ_SRC_DIRS    - Put a list of sub-directories with code that should be 
#                   compiled into explicit .o files but NOT stored in the 
#                   archive library here (for instance, object files for 
#                   for programs).
# CONFIG_DIRS     - Put a list of sub-directories with initialization files
#                   here.   NOTE:  It is assumed that ALL config files in the
#                   CONFIG_DIRS are of the form *.*  This is to avoid copying
#                   of contained sub-directories (eg, the CVS sub-directory).
# LOCAL_INCS      - Local directories to search for include files
# LOCAL_MODS      - Local directories to search for F90 compiled module files
#-----------------------------------------------------------------------------
LIB_SRC_DIRS := ./code ./custom ./cesr ./modules ./cpp_interface
OBJ_SRC_DIRS := ./test_jobs

CONFIG_DIRS   :=
LOCAL_INCS    := ./include ../include
LOCAL_MODS    := ../modules 

#SRC_EXCLUDE := ele_to_fibre.f90

#-----------------------------------------------------------------------------
# "EXTRA" variables can be specified at the command line or hardwired here.  
# These variables are automatically appended to the relevant search paths and
# lists in M.tail
#
# EXTRA_SRC      - extra code source directories
# EXTRA_MOD_SRC  - extra module source directories
# EXTRA_OBJ_SRC  - extra object source directories
# EXTRA_CONFIG   - extra configuration file directories
# EXTRA_INCS     - extra include file search directories
# EXTRA_MODS     - extra compile module search directories
# EXTRA_FFLAGS   - extra Fortran flags for compilation
# EXTRA_CFLAGS   - extra C flags for compilation 
# EXTRA_LFLAGS   - extra linker flags
#-----------------------------------------------------------------------------

#------------------------------------------------
# Name of local libraries (standard and debug) 
#------------------------------------------------
ifeq "$(LIBNAME)" ""
  WHERE   := $(shell pwd)
  LIBNAME := $(notdir $(WHERE))
endif


#------------------------------------------------
# Include the standard CESR M.tail makefile
#------------------------------------------------
ifeq "$(BMAD_DIST)" "TRUE"
  include $(BMAD_GMAKE)/M.tail
else
  include $(CESR_GMAKE)/M.tail
endif


