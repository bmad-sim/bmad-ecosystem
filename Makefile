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
# $Id$
#
# $Log$
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

#------------------------------------------------
# Set up some standard locally referenced paths
#------------------------------------------------
absrcdir = $(shell pwd)


#---------------------------------------------------------------------------
# Specify default source directories for code to be compiled and put into 
# libraries, for module files, for code to be compiled into explicit object
# files, and initialization files to be copied to the config directory tree.  
#
# Multiple directories should be separated by spaces.
#
# DO NOT put a trailing SLASH on directory names
#---------------------------------------------------------------------------
CODE_SRC_DIRS := ./code
MOD_SRC_DIRS  := ./modules
OBJ_SRC_DIRS  := 
# It is assumed that ALL config files in the CONFIG_DIRS are of the form *.*.  
# This is to avoid copying of directories (eg, the CVS sub-directory) 
CONFIG_DIRS   := ./lat ./layout


#------------------------------------------------
# Specify directories to be searched for include
# files and compiled f90 module files.
#------------------------------------------------
LOCAL_INCS := ./include ../include
OTHER_INCS := $(CESR_INC) $(EXTRA_INCS)

LOCAL_MODS := ../modules
OTHER_MODS := $(CESR_MOD) $(PKG_MOD) $(EXTRA_MODS)


#------------------------------------------------
# Specify output directories for archive 
# libraries and compiled module files
#------------------------------------------------
locallib    := ../lib
localmod    := ../modules
localconfig := ../config


#------------------------------------------------
# Name of local libraries (standard and debug) 
#------------------------------------------------
ifeq "$(LIBNAME)" ""
  WHERE   := $(shell pwd)
  LIBNAME := $(notdir $(WHERE))
endif


#------------------------------------------------
# Set up production versus debug compilation
#------------------------------------------------
ifneq (,$(findstring $(DEBUG),yY))
  LIBRARY     := $(locallib)/lib$(LIBNAME)_g.a
  MOD_OUT_DIR := $(localmod)_g
  FFLAGS      := -g
  CCFLAGS     := -g
else
  LIBRARY       := $(locallib)/lib$(LIBNAME).a
  MOD_OUT_DIR   := $(localmod)
  LIBRARY_G     := $(locallib)/lib$(LIBNAME)_g.a
  MOD_OUT_DIR_G := $(localmod)_g
  FFLAGS        :=
  CCFLAGS       :=
endif


#------------------------------------------------
# Set up some standard definitions
#------------------------------------------------
SHELL := /bin/sh
FC := f90 -cpp
CC := /usr/bin/cc
DEBUG := n
MAKE := gmake
LIB_MAKEFILE := ./Makefile
CP := /usr/local/bin/cp
RM := /usr/local/bin/rm -f
MV := /usr/local/bin/mv -v
CD := /usr/bin/cd
# special cesr fortran-dependency script
CESRDEP := $(CESR_UTIL)/cesrdepend


#------------------------------------------------
# Set up some standard compilation and linking 
# flags
#------------------------------------------------
FFLAGS += -c -u -v -module $(MOD_OUT_DIR)
CFLAGS += -c -check
LFLAGS := -Wl,-m -non_shared


#------------------------------------------------
# Generate the full CESR include file and module 
# search path.
#------------------------------------------------
LIST_INCS  := $(patsubst %, -I%, $(LOCAL_INCS))
LIST_INCS  += $(patsubst %, -I%, $(OTHER_INCS))
LIST_INCS  += $(patsubst %, -I%, $(LOCAL_MODS))
LIST_INCS  += $(patsubst %, -I%, $(OTHER_MODS))


#------------------------------------------------
# Dependency file
#------------------------------------------------
DEP_FILE := $(absrcdir)/.depend


#---------------------------------------------------------------------------
# Generate list of fortran and c files for compilation.  
# Fortran sources may have the following suffixes:  .f90, .f77, .for, and .F
#---------------------------------------------------------------------------
F90_SEARCH             := $(patsubst %, %/*.f90, $(CODE_SRC_DIRS))
F90_SOURCES            := $(wildcard $(F90_SEARCH))
NOTDIR_F90_SOURCES     := $(notdir $(F90_SOURCES))
NOTDIR_F90_OBJS        := $(NOTDIR_F90_SOURCES:%.f90=$(LIBRARY)(%.o))

F77_SEARCH             := $(patsubst %, %/*.f77, $(CODE_SRC_DIRS))
F77_SOURCES            := $(wildcard $(F77_SEARCH))
NOTDIR_F77_SOURCES     := $(notdir $(F77_SOURCES))
NOTDIR_F77_OBJS        := $(NOTDIR_F77_SOURCES:%.f77=$(LIBRARY)(%.o))

FOR_SEARCH             := $(patsubst %, %/*.for, $(CODE_SRC_DIRS))
FOR_SOURCES            := $(wildcard $(FOR_SEARCH))
NOTDIR_FOR_SOURCES     := $(notdir $(FOR_SOURCES))
NOTDIR_FOR_OBJS        := $(NOTDIR_FOR_SOURCES:%.for=$(LIBRARY)(%.o))

F_SEARCH               := $(patsubst %, %/*.F, $(CODE_SRC_DIRS))
F_SOURCES              := $(wildcard $(F_SEARCH))
NOTDIR_F_SOURCES       := $(notdir $(F_SOURCES))
NOTDIR_F_OBJS          := $(NOTDIR_F_SOURCES:%.F=$(LIBRARY)(%.o))

FORTRAN_SOURCES        := $(F90_SOURCES) 
FORTRAN_SOURCES        += $(F77_SOURCES) 
FORTRAN_SOURCES        += $(FOR_SOURCES) 
FORTRAN_SOURCES        += $(F_SOURCES)

NOTDIR_FORTRAN_SOURCES := $(NOTDIR_F90_SOURCES) 
NOTDIR_FORTRAN_SOURCES += $(NOTDIR_F77_SOURCES) 
NOTDIR_FORTRAN_SOURCES += $(NOTDIR_FOR_SOURCES) 
NOTDIR_FORTRAN_SOURCES += $(NOTDIR_F_SOURCES)

NOTDIR_FORTRAN_OBJS    := $(NOTDIR_F90_OBJS) 
NOTDIR_FORTRAN_OBJS    += $(NOTDIR_F77_OBJS) 
NOTDIR_FORTRAN_OBJS    += $(NOTDIR_FOR_OBJS) 
NOTDIR_FORTRAN_OBJS    += $(NOTDIR_F_OBJS)

C_SEARCH               := $(patsubst %, %/*.c, $(CODE_SRC_DIRS))
C_SOURCES              := $(wildcard $(C_SEARCH))
NOTDIR_C_SOURCES       := $(notdir $(C_SOURCES))
NOTDIR_C_OBJS          := $(NOTDIR_C_SOURCES:%.c=$(LIBRARY)(%.o))

NOTDIR_ALL_OBJS        := $(NOTDIR_FORTRAN_OBJS) $(NOTDIR_C_OBJS)


#------------------------------------------------
# Generate list of module sources for compilation
#------------------------------------------------
MODULE_SEARCH      := $(patsubst %, %/*.f90, $(MOD_SRC_DIRS))
MODULE_SOURCES     := $(wildcard $(MODULE_SEARCH))
NOTDIR_MOD_SOURCES := $(notdir $(MODULE_SOURCES))
NOTDIR_MODULES     := $(NOTDIR_MOD_SOURCES:%.f90=%.mod)
NOTDIR_MOD_OBJS    := $(NOTDIR_MOD_SOURCES:%.f90=$(LIBRARY)(%.o))
OUTPUT_MODULES     := $(NOTDIR_MODULES:%.mod=$(MOD_OUT_DIR)/%.mod)


#------------------------------------------------
# Set up full source and dependency lists
#------------------------------------------------
ALL_SOURCES := $(FORTRAN_SOURCES) $(C_SOURCES) $(MODULE_SOURCES) 
ALL_DEPEND  := $(F90_SOURCES:.f90=.d) 
ALL_DEPEND  += $(F77_SOURCES:.f77=.d) 
ALL_DEPEND  += $(FOR_SOURCES:.for=.d) 
ALL_DEPEND  += $(F_SOURCES:.F=.d) 
ALL_DEPEND  += $(C_SOURCES:.c=.d) 
ALL_DEPEND  += $(MODULE_SOURCES:.f90=.dmd)


#----------------------------------------------------------
# Generate the configuration sources and final destinations
#----------------------------------------------------------
CONFIG_OUT_DIR  := $(localconfig)/$(LIBNAME)
CONFIG_OUT_SUBS := $(patsubst %, $(CONFIG_OUT_DIR)/%, $(notdir $(CONFIG_DIRS)))

CONFIG_SEARCH  := $(patsubst %, %/*.*, $(CONFIG_DIRS))
CONFIG_SOURCES := $(wildcard $(CONFIG_SEARCH))
CONFIG_OUTPUT  := $(patsubst %, $(CONFIG_OUT_DIR)/%, $(CONFIG_SOURCES))


#------------------------------------------------------------
# here are some sneaky make tricks to unbuild the UNDEF list
# (Thanks to M Athanas)
#------------------------------------------------------------
COMMA := ,
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)

VPATH_LIST := . $(CODE_SRC_DIRS) $(MOD_SRC_DIRS)
VPATH      := $(subst $(SPACE),:,$(VPATH_LIST))


#------------------------------------------------
# Set up actions
#------------------------------------------------
.PHONY: all modules_only mklib mkmod mkconfig production debug modules modules_debug config clean dump depend
	@

all: depend mklib mkmod mkconfig production debug config

modules_only: depend mklib mkmod modules modules_debug

mklib:
	@test -z "$(LIBRARY)" || test -d $(locallib) || mkdir -p $(locallib)
	@test -z "$(LIBRARY)" || test -f $(LIBRARY)  || ar r $(LIBRARY) ;
	@test -z "$(LIBRARY_G)" || test -f $(LIBRARY_G)  || ar r $(LIBRARY_G)

mkmod:
	@test -d $(MOD_OUT_DIR) || mkdir -p $(MOD_OUT_DIR)
	@test -d $(MOD_OUT_DIR_G) || mkdir -p $(MOD_OUT_DIR_G)

mkconfig:
	@test -z "$(CONFIG_DIRS)" || test -d $(localconfig) || mkdir -p $(localconfig)
	@test -z "$(CONFIG_DIRS)" || test -d $(CONFIG_OUT_DIR)  || mkdir -p $(CONFIG_OUT_SUBS)

production: modules $(LIBRARY)

debug: 
	@$(MAKE) DEBUG=y production

modules_debug: 
	@$(MAKE) DEBUG=y modules

modules: $(OUTPUT_MODULES)

# Make library dependent on output modules since the object files 
# are installed in the library when the modules are built
$(LIBRARY): $(NOTDIR_ALL_OBJS) $(OUTPUT_MODULES)
	@echo " "
	@echo "Object file table for:  $(LIBRARY)"
	@$(AR) ts $(LIBRARY)
	@echo " "

config: $(CONFIG_OUTPUT)

clean:
	@$(RM) $(LIBRARY)
	@$(RM) $(LIBRARY:.a=_g.a)
	@$(RM) $(NOTDIR_MODULES:%.mod=$(MOD_OUT_DIR)/%.mod)
	@$(RM) $(NOTDIR_MODULES:%.mod=$(MOD_OUT_DIR_G)/%.mod)
	@$(RM) $(DEP_FILE)
	@$(RM) $(ALL_DEPEND)

dump:
	@echo "STD Library is $(LIBRARY)"
	@echo "DBG Library is $(LIBRARY_G)"
	@echo "Compiled module output to $(MOD_OUT_DIR)"
	@echo "Debug    module output to $(MOD_OUT_DIR_G)"
	@echo "Code source directories are:  "
	@echo " $(CODE_SRC_DIRS)"
	@echo "Fortran sources are: "
	@echo " $(NOTDIR_FORTRAN_SOURCES)"
	@echo "C sources are: "
	@echo " $(C_SOURCES)"
	@echo "Module source directories are:"
	@echo " $(MOD_SRC_DIRS)"
	@echo "Modules are: "
	@echo " $(NOTDIR_MODULES)"  
	@echo "Configuration file directories are:"
	@echo " $
	@echo "Include path list is:"
	@echo " $(LIST_INCS)"	
	@echo "Search path list is:"
	@echo " $(VPATH)"	

depend: $(DEP_FILE)
	@

.INTERMEDIATE: $(ALL_DEPEND)
	@

$(DEP_FILE):  $(ALL_DEPEND)
	@test -z "$(DEP_FILE)" || test -f $(DEP_FILE)  || touch $(DEP_FILE)


#------------------------------------------------
# Set up pattern rules
#------------------------------------------------
%.d : %.f90
	@test -f $@ || touch $@
# Always put in the library target
	@cd $(dir $<); makedepend $(LIST_INCS) -o".o) $(notdir $@)" -p$$\(LIBRARY\)\( -f$(notdir $@) $(notdir $<)
# Look for fortran module and include dependencies
	@$(CESRDEP) $(LIST_INCS) -MOD"$(NOTDIR_MODULES)" -MP"$$(MOD_OUT_DIR)" -SUF".o) $(notdir $@)" -PRE"$$(LIBRARY)(" $< >> $@
	@cat $@ | grep -v "DO NOT DELETE" >> $(DEP_FILE)		
	@$(RM) $@.bak

%.d : %.f77
	@test -f $@ || touch $@
# Always put in the library target
	@cd $(dir $<); makedepend $(LIST_INCS) -o".o) $(notdir $@)" -p$$\(LIBRARY\)\( -f$(notdir $@) $(notdir $<)
# Look for fortran module and include dependencies
	@$(CESRDEP) $(LIST_INCS) -MOD"$(NOTDIR_MODULES)" -MP"$$(MOD_OUT_DIR)" -SUF".o) $(notdir $@)" -PRE"$$(LIBRARY)(" $< >> $@
	@cat $@ | grep -v "DO NOT DELETE" >> $(DEP_FILE)		
	@$(RM) $@.bak

%.d : %.for
	@test -f $@ || touch $@
# Always put in the library target
	@cd $(dir $<); makedepend $(LIST_INCS) -o".o) $(notdir $@)" -p$$\(LIBRARY\)\( -f$(notdir $@) $(notdir $<)
# Look for fortran module and include dependencies
	@$(CESRDEP) $(LIST_INCS) -MOD"$(NOTDIR_MODULES)" -MP"$$(MOD_OUT_DIR)" -SUF".o) $(notdir $@)" -PRE"$$(LIBRARY)(" $< >> $@
	@cat $@ | grep -v "DO NOT DELETE" >> $(DEP_FILE)		
	@$(RM) $@.bak

%.d : %.F
	@test -f $@ || touch $@
# Always put in the library target
	@cd $(dir $<); makedepend $(LIST_INCS) -o".o) $(notdir $@)" -p$$\(LIBRARY\)\( -f$(notdir $@) $(notdir $<)
# Look for fortran module and include dependencies
	@$(CESRDEP) $(LIST_INCS) -MOD"$(NOTDIR_MODULES)" -MP"$$(MOD_OUT_DIR)" -SUF".o) $(notdir $@)" -PRE"$$(LIBRARY)(" $< >> $@
	@cat $@ | grep -v "DO NOT DELETE" >> $(DEP_FILE)		
	@$(RM) $@.bak

%.dmd : %.f90
	@test -f $@ || touch $@
# Produce dependencies for both the output .mod file as well as the library object file
	@cd $(dir $<); makedepend $(LIST_INCS) -o".mod $(notdir $@)" -p$$\(MOD_OUT_DIR\)\/ -f$(notdir $@) $(notdir $<)
	@$(CESRDEP) $(LIST_INCS) -MOD"$(NOTDIR_MODULES)" -MP"$$(MOD_OUT_DIR)" -SUF".mod $(notdir $@)" -PRE"$$(MOD_OUT_DIR)/" $< >> $@
	@cd $(dir $<); makedepend -a $(LIST_INCS) -o".o) $(notdir $@)" -p$$\(LIBRARY\)\( -f$(notdir $@) $(notdir $<)
	@$(CESRDEP) $(LIST_INCS) -MOD"$(NOTDIR_MODULES)" -MP"$$(MOD_OUT_DIR)" -SUF".o) $(notdir $@)" -PRE"$$(LIBRARY)(" $< >> $@
	@cat $@ | grep -v "DO NOT DELETE" >> $(DEP_FILE)		
	@$(RM) $@.bak

%.d : %.c
	@test -f $@ || touch $@
	@cd $(dir $<); makedepend $(LIST_INCS) -o".o) $(notdir $@)" -p$$\(LIBRARY\)\( -f$(notdir $@) $(notdir $<)
	@cat $@ | grep -v "DO NOT DELETE" >> $(DEP_FILE)		
	@$(RM) $@.bak

(%.o) : %.f90
#	@echo "Compiling:  $@ $%"
	$(FC)  $(FFLAGS) $(LIST_INCS) $<
	$(AR) r $@ $%
	@$(RM) $%

(%.o) : %.f77
#	@echo "Compiling:  $@ $%"
	$(FC)  $(FFLAGS) $(LIST_INCS) $<
	$(AR) r $@ $%
	@$(RM) $%

(%.o) : %.for
#	@echo "Compiling:  $@ $%"
	$(FC)  $(FFLAGS) $(LIST_INCS) $<
	$(AR) r $@ $%
	@$(RM) $%

(%.o) : %.F
#	@echo "Compiling:  $@ $%"
	$(FC)  $(FFLAGS) $(LIST_INCS) $<
	$(AR) r $@ $%
	@$(RM) $%

(%.o) : %.c
#	@echo "Compiling:  $@ $%"
	$(CC)  $(CFLAGS) $(LIST_INCS) $<
	$(AR) r $@ $%
	@$(RM) $%

$(MOD_OUT_DIR)/%.mod : %.f90
#	@echo "Compiling:  $@ $<"
	$(FC)  $(FFLAGS) $(LIST_INCS) $<
	$(AR) r $(LIBRARY) $(notdir $(@:.mod=.o))
	@$(RM) $(notdir $(@:.mod=.o))

$(CONFIG_OUT_DIR)/% : %
	$(CP) $< $(CONFIG_OUT_DIR)/$<

%.inc : 
	@


#------------------------------------------------
#  Get rid of some useless rules
#------------------------------------------------
%::%,v
%::RCS/%,v
%::s.%
%::SCCS/s.%
%::%.S
%::%.s
%::%.web
%::%.p
%::%.l


#------------------------------------------------
# Include dependency file (- option means don't 
# abort if not available)
#------------------------------------------------
-include $(DEP_FILE)



