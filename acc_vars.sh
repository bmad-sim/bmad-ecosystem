#--------------------------------------------------------------
# acc_vars.sh
#
# ACCelerator libraries management and software build system
# environment setup script.
#
# This script uses Bourne shell syntax. When sourced into a
# user's environment it allows for selection of a particular
# pre-built accelerator libraries release within the lab's
# Linux computing landscape.
#
# Variables that can be used to control the behavior of this
# environment setup process:
#   To set the release to use when searching for include
#   files and when linking code, several variable names are
#   honored here if they exist in the sourcing user's 
#   environment to accommodate historical use habits and old
#   lab-wide conventions.
#
#    CESRLIB is honored unless
#    ACCLIB is set, which is honored unless
#    ACC_RELEASE_REQUEST is set.
#
#  Environment setup is configured, for the time-being, to
#  force 32-bit building and linking on all systems,
#  regardless of machine type/architecture.  To change this
#  default behavior, the user must set the variable
#
#    ACC_FORCE_32_BIT    "N" 
#
#  in their environment BEFORE sourcing the setup script.
#  This variable is only honored if the machine on which
#  the setup script is sourced allows the option of running
#  both 32-bit and 64-bit code.  Otherwise, on a 32-bit-only
#  system it is simply ignored.
#
#
# Once this script is sourced, ideally by getting called from
# the user's login script, shortcut aliases are available to
# allow easy selection of an active release and architecture.
#  'current' -- Sets active relase to CURRENT
#  'devel' -- Sets active release to DEVEL
#  'setrel <release_name> -- Sets a specific release as active
#
# Other useful aliases defined below are:
#  'accinfo' -- Displays a summary of the active release's
#      short name ('devel', 'current' etc), it's true (dated)
#      name, build architecture, and selected Fortran compiler.
#  
# For C-shell variant use, all of the above applies since the
# C-shell script is a wrapper for this one.  The user must
# source the script named 
#   acc_vars.csh  instead.
#--------------------------------------------------------------

#--------------------------------------------------------------
# For User defined compiler, check how ACC_SET_F_COMPILER
# is defined, if not defined then set to default "ifort"
#--------------------------------------------------------------
if ( [ -z ${ACC_SET_F_COMPILER} ] || [ "${ACC_SET_F_COMPILER}" != "gfortran" ]  ) then
    export ACC_SET_F_COMPILER="ifort"
fi

#--------------------------------------------------------------

ONLINE_ARCHIVE_BASE_DIR='/nfs/cesr/online/lib'
ONLINE_OPT_DIR='/nfs/cesr/opt'

ONLINE_RELEASE_MGMT_DIR=${ONLINE_ARCHIVE_BASE_DIR}'/util'

#ONLINE_IFORT_SETUP_DIR=${ONLINE_OPT_DIR}'/intel/parallel_studio_xe_2016.2.062/bin'
ONLINE_IFORT_SETUP_DIR=${ONLINE_OPT_DIR}'/intel/parallel_studio_xe_2017_update7/bin'
ONLINE_IFORT_SETUP_COMMAND=(${ONLINE_IFORT_SETUP_DIR}'/compilervars.sh intel64')

ONLINE_GFORTRAN_SETUP_DIR='/opt/rh/devtoolset-4/'
ONLINE_GFORTRAN_SETUP_COMMAND=(${ONLINE_GFORTRAN_SETUP_DIR}/enable)

#--------------------------------------------------------------

OFFLINE_ARCHIVE_BASE_DIR='/nfs/acc/libs'
OFFLINE_OPT_DIR='/nfs/opt'

[ "${OFFLINE_LOCAL_ARCHIVE_BASE_DIR}" ] && OFFLINE_ARCHIVE_BASE_DIR=${OFFLINE_LOCAL_ARCHIVE_BASE_DIR}

OFFLINE_RELEASE_MGMT_DIR=${OFFLINE_ARCHIVE_BASE_DIR}'/util'

#OFFLINE_IFORT_SETUP_DIR=${OFFLINE_OPT_DIR}'/intel/parallel_studio_xe_2016.2.062/bin'
OFFLINE_IFORT_SETUP_DIR=${OFFLINE_OPT_DIR}'/intel/parallel_studio_xe_2017_update7/bin'
OFFLINE_IFORT_SETUP_COMMAND=(${OFFLINE_IFORT_SETUP_DIR}'/compilervars.sh intel64')

OFFLINE_GFORTRAN_SETUP_DIR='/opt/rh/devtoolset-4/'
OFFLINE_GFORTRAN_SETUP_COMMAND=(${OFFLINE_GFORTRAN_SETUP_DIR}/enable)

# Capture value of ACC_BIN to allow removal from path for cleanliness.
OLD_ACC_BIN=${ACC_BIN}

#--------------------------------------------------------------
# Support functions
#--------------------------------------------------------------
# Is argument $1 missing from argument $2 (or PATH)?
function no_path() {
    eval "case :\$${2-PATH}: in *:$1:*) return 1;; *) return 0;; esac"
}

# If argument path ($1) exists and is not in path, append it.
function add_path() {
    [ -d ${1:-.} ] && no_path $* && eval ${2:-PATH}="\$${2:-PATH}:$1"
}

# If argument ($1) is in PATH (or 2nd argument), remove it.
function del_path() {
    no_path $* || eval ${2:-PATH}=`eval echo :'$'${2:-PATH}: |
    /bin/sed -e "s;:$1:;:;g" -e "s;^:;;" -e "s;:\$;;"`
}

# Test whether argument 1 string contains argument 2 as a substring.
# Returns true(0) if substring is present
# Returns false(1) is substring is not
function has_substring() {
    if [[ "${1}" =~ "${2}" ]]; then
	return `true`
    else
	return `false`
    fi
}

# Function to remove duplicate path entries for any path variable
# passed as an argument.
function remove_duplicates() {
    echo ${1} | awk -F: '{for (i=1;i<=NF;i++) { if ( !x[$i]++ ) printf("%s:",$i); }}'
}


#--------------------------------------------------------------
# Files to hold a 'printenv' environment snapshot before this
# script defines all its variables, and again after for use
# in the wrapper script that only C-style shells will source.
# This is so that all the syntax for environment setup remains
# in one language, (bash) and need not be kept up-to-date
# between a bourne-shell setup script and C-shell one.
#
# Take first environment snapshot for C-shell support.
#--------------------------------------------------------------
if ( [ "${ENV_USE_SNAPSHOTS}" == "Y" ] ) then
    ENV_SNAPSHOT_1="${HOME}/.ACC_env_snapshot_1.tmp"
    ENV_SNAPSHOT_2="${HOME}/.ACC_env_snapshot_2.tmp"
    printenv > ${ENV_SNAPSHOT_1}
fi


#--------------------------------------------------------------
# Set CESR_ONLINE to the path used for 'CESR offline'
# hosts if it hasn't already been set upon entering this
# script.
#--------------------------------------------------------------
export CESR_ONLINE=${CESR_ONLINE:-/nfs/cesr/online}



#--------------------------------------------------------------
# 'CESR ONLINE' libraries are kept up-to-date with /nfs/acc/libs
# in case of a network failure.  In that event, change the 
# following variable to 'true' to use the online libs.
#--------------------------------------------------------------

case $(uname -n) in
    cesr*) USE_CESR_ONLINE_LOC='true';;
        *) USE_CESR_ONLINE_LOC='false';;
esac



if ( [ "${USE_CESR_ONLINE_LOC}" == "true" ] ) then

    #-- CESR Online Host
    SETUP_SCRIPTS_DIR=${ONLINE_RELEASE_MGMT_DIR}
    RELEASE_ARCHIVE_BASE_DIR=${ONLINE_ARCHIVE_BASE_DIR}
    HOST_TYPE="CESR ONLINE"

    

else
    #-- CESR OFFLINE Host
    SETUP_SCRIPTS_DIR=${OFFLINE_RELEASE_MGMT_DIR}
    RELEASE_ARCHIVE_BASE_DIR=${OFFLINE_ARCHIVE_BASE_DIR}
    HOST_TYPE="CESR OFFLINE"

fi 

# Online computers should use the online compiler and debugger
case $(uname -n) in
    cesr*) IFORT_SETUP_COMMAND=${ONLINE_IFORT_SETUP_COMMAND}
	GFORTRAN_SETUP_COMMAND=${ONLINE_GFORTRAN_SETUP_COMMAND}
	export LM_LICENSE_FILE="/opt/intel/licenses/l_N6CKW3X8.lic:${ONLINE_OPT_DIR}/totalview/license.dat:${OFFLINE_OPT_DIR}/totalview/license.dat:${ONLINE_OPT_DIR}/MathWorks/Matlab_r2013a/licenses/network.lic:${OFFLINE_OPT_DIR}/MathWorks/Matlab_r2013a/licenses/network.lic"
	export PATH=${PATH}:${ONLINE_OPT_DIR}/totalview/bin
	;;					 
    *) IFORT_SETUP_COMMAND=${OFFLINE_IFORT_SETUP_COMMAND}
        GFORTRAN_SETUP_COMMAND=${OFFLINE_GFORTRAN_SETUP_COMMAND}
	export LM_LICENSE_FILE="${OFFLINE_OPT_DIR}/totalview/license.dat:${ONLINE_OPT_DIR}/totalview/license.dat:${OFFLINE_OPT_DIR}/MathWorks/Matlab_r2013a/licenses/network.lic:${ONLINE_OPT_DIR}/MathWorks/Matlab_r2013a/licenses/network.lic"
	export PATH=${PATH}:${OFFLINE_OPT_DIR}/totalview/bin
	;;					 
esac




#--------------------------------------------------------------
# Allow for util directory override.
#--------------------------------------------------------------
SETUP_SCRIPTS_DIR=${UTIL_DIR_REQUEST:-$SETUP_SCRIPTS_DIR}
unset UTIL_DIR_REQUEST


#--------------------------------------------------------------
# Variables that influence the build process 
# (Their names begin with 'ACC'.)
#--------------------------------------------------------------
ACC_OS="`uname`"
export ACC_ARCH="`uname -m`"
export ACC_OS_ARCH="${ACC_OS}_${ACC_ARCH}"

if ( [ "${ACC_SET_F_COMPILER}" == "ifort" ] ) then
    export ACC_FC="intel"
else
    export ACC_FC="gfortran"
fi


#--------------------------------------------------------------
# Allow for default behavior on 64-bit systems to be
# overridden to produce and use 32-bit binaries.  
#
# Until full 64-bit OS migration takes place, all builds will
# default to 32-bit on all hosts.  This can be changed to 
# use the true machine architecture word length per user
# request by setting the ACC_FORCE_32_BIT variable to the 
# value "N".
#--------------------------------------------------------------
if ( [ "${ACC_FORCE_32_BIT}" == "" ] || [ "${ACC_FORCE_32_BIT}" == "N" ] ) then
    #ACC_FORCE_32_BIT=N
    export ACC_ARCH="x86_64"
    export ACC_OS_ARCH="${ACC_OS}_${ACC_ARCH}"
fi

if (  [ "${ACC_FORCE_32_BIT}" = "Y" ] ) then
    export ACC_ARCH="i686"
    export ACC_OS_ARCH="${ACC_OS}_${ACC_ARCH}"
fi



export ACC_PLATFORM="${ACC_OS_ARCH}_${ACC_FC}"
PLATFORM_DIR=${RELEASE_ARCHIVE_BASE_DIR}/${ACC_PLATFORM}

export ACC_BASE=${PLATFORM_DIR}



#--------------------------------------------------------------
# Determine if the requested release name was specified, AND
# if the that subdirectory of the platform directory whose 
# name was generated above exists.
# If not, use the default "stable" release.
#
# Several variable names are honored here if they exist in the
# sourcing user's environment to accommodate historical use
# habits.
#    CESRLIB is honored unless
#    ACCLIB is set, which is honored unless
#    ACC_RELEASE_REQUEST is set. 
#--------------------------------------------------------------
if ( [ "${CESRLIB}" != "" ] ) then
    RELEASE_REQUEST=${CESRLIB}
fi
if ( [ "${ACCLIB}" != "" ] ) then
    RELEASE_REQUEST=${ACCLIB}
fi
if ( [ "${ACC_RELEASE_REQUEST}" != "" ] ) then
    RELEASE_REQUEST=${ACC_RELEASE_REQUEST}
fi

if ( [ "${RELEASE_REQUEST}" != "" ] ) then
    if ( [ "${RELEASE_REQUEST}" != "" -a -d ${PLATFORM_DIR}/${RELEASE_REQUEST} ] ) then
        export ACC_RELEASE=${RELEASE_REQUEST}
    else
        echo "Release request \"${RELEASE_REQUEST}\" not found in"
        echo "${PLATFORM_DIR}.  Defaulting to \"current\" release."
        export ACC_RELEASE="current"
    fi
else
    export ACC_RELEASE="current"
fi
unset CESRLIB
unset ACCLIB
unset ACC_RELEASE_REQUEST

export ACC_RELEASE_DIR=${PLATFORM_DIR}/${ACC_RELEASE}
export ACC_ROOT_DIR=${ACC_RELEASE_DIR} # Defined for both Release and Distribution Builds
export ACC_TRUE_RELEASE=`readlink ${ACC_RELEASE_DIR}`
if ( [ "${ACC_TRUE_RELEASE}" == "" ] ) then
    ACC_TRUE_RELEASE=${ACC_RELEASE}
fi
export ACC_UTIL=${SETUP_SCRIPTS_DIR}
export ACC_BIN=${ACC_RELEASE_DIR}/production/bin
export ACC_EXE=${ACC_BIN} # For backwards compatibility
export ACC_DEBUG=${ACC_RELEASE_DIR}/debug/bin
export ACC_PKG=${ACC_RELEASE_DIR}/packages
export ACC_REPO=https://accserv.lepp.cornell.edu/svn/
export ACCR=https://accserv.lepp.cornell.edu/svn/
export ACC_BUILD_SYSTEM=${RELEASE_ARCHIVE_BASE_DIR}/build_system
export ACC_BUILD_EXES=Y
export ACC_CMAKE_VERSION=3.13.2 # Updated as per request in RT#56874
export TAO_DIR=${ACC_RELEASE_DIR}/tao

#--------------------------------------------------------------
# Variables to support other software packages
#--------------------------------------------------------------
export PGPLOT_DIR=${ACC_PKG}/PGPLOT
export PGPLOT_FONT=${ACC_PKG}/PGPLOT/grfont.dat



#--------------------------------------------------------------
# Set up the fortran compiler appropriate to the machine
# architecture.
#
# TODO: These vendor-supplied compiler setup scripts will
#       need to live in both CESR OFFLINE and CESR ONLINE 
#       directories.
#--------------------------------------------------------------
case ${ACC_OS_ARCH} in

    "Linux_x86_64" )
	if ( [ "${ACC_SET_F_COMPILER}" == "ifort" ] ) then
	    source ${IFORT_SETUP_COMMAND}
	else
	    source ${GFORTRAN_SETUP_COMMAND}
	fi
	;;

esac



#--------------------------------------------------------------
# Append necessary directories to the user's PATH environment
# variable.
#--------------------------------------------------------------
del_path ${OLD_ACC_BIN}

add_path ${ACC_UTIL}
add_path ${ACC_PKG}/bin
add_path ${ACC_BIN}


#--------------------------------------------------------------
# Append path to guarantee proper cmake executable is used
# Added backup PATH /nfs/opt/extra/bin as requested in RT#65470
#--------------------------------------------------------------

if [ -d ${PLATFORM_DIR}/extra/bin ]; then
	PATH=${PLATFORM_DIR}/extra/bin:$PATH
elif [ -d /nfs/opt/extra/bin ]; then
	PATH=/nfs/opt/extra/bin:$PATH
fi

#--------------------------------------------------------------
# Ensure that there is only one reference to the release's
# lib directory by searching the LD_LIBRARY_PATH variable
# value for paths that are prefixed with the
# RELEASE_ARCHIVE_BASE_DIR and removing them.  The last
# step is appending the presently valid lib path for
# the active release.
#--------------------------------------------------------------
DEFAULT_IFS=${IFS} # bash 'internal field separator' reserved variable
IFS=":"
LD_LIBRARY_PATH_TEMP=""
for part in ${LD_LIBRARY_PATH}; do
  case ${part} in

      "${RELEASE_ARCHIVE_BASE_DIR}"*)
	#  echo "ACC path found. Will remove from LD_LIBRARY_PATH."
	  continue
	  ;;

      #"")
      #    continue
#	  ;;

      *)
        if ( [ "${LD_LIBRARY_PATH_TEMP}" == "" ] ) then
	  LD_LIBRARY_PATH_TEMP=${part}
	else
	  LD_LIBRARY_PATH_TEMP=${LD_LIBRARY_PATH_TEMP}:${part}
	fi
	;;

  esac
done
IFS=${DEFAULT_IFS}
LD_LIBRARY_PATH=${LD_LIBRARY_PATH_TEMP}


#--------------------------------------------------------------
# Source ROOT setup script for the set of ROOT libraries and
# programs built as 3rd-party packages in the selected release.
#--------------------------------------------------------------
#if [ -f ${ACC_PKG}/root/bin/thisroot.sh ]; then
    # To prevent multiple sourcing
    #has_substring ${PATH} "root/../bin" 
    #if [ "${?}" == "1" ]; then
	source ${ACC_PKG}/production/bin/thisroot.sh
    #fi
#fi

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ACC_RELEASE_DIR}/production/lib:${ACC_RELEASE_DIR}/packages/production/lib


#--------------------------------------------------------------
# Useful aliases
#--------------------------------------------------------------
# Information about build system setup for the current shell session.
function accinfo() {
    printenv | grep ACC_ | sort;
    export ACC_TRUE_RELEASE=`readlink ${ACC_RELEASE_DIR}`
    echo -e "ACC release      : \"${ACC_RELEASE}\" [${ACC_TRUE_RELEASE}]
Architecture     : ${ACC_ARCH}
Fortran compiler : ${ACC_FC}
Cmake version    : $(cmake --version | head -1 | awk ' { print $3 } ')"
}

alias ACCINFO='accinfo'
alias acc_info='accinfo'
alias ACC_INFO='accinfo'


alias current='ACC_RELEASE_REQUEST=current; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'
alias devel='ACC_RELEASE_REQUEST=devel; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'
alias nightly='ACC_RELEASE_REQUEST=nightly; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'
#alias online='ACC_RELEASE_REQUEST=online; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'
#alias online-devel='ACC_RELEASE_REQUEST=online-devel; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'

alias current32='ACC_FORCE_32_BIT=Y; current'
alias devel32='ACC_FORCE_32_BIT=Y; devel;'

# These are functions so that they may take arguments from the user.

function setrel() {
    ACC_RELEASE_REQUEST=${1}; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh
}

function setrel64() {
    ACC_FORCE_32_BIT=N; ACC_RELEASE_REQUEST=${1}; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh
}


#--------------------------------------------------------------
# Remove all duplicate path entries from variables composed during
# the above process
#--------------------------------------------------------------
NEWPATH=$(remove_duplicates ${PATH})
export PATH=${NEWPATH%?}
export LD_LIBRARY_PATH=`remove_duplicates ${LD_LIBRARY_PATH}`
export DYLD_LIBRARY_PATH=`remove_duplicates ${DYLD_LIBRARY_PATH}`
export MANPATH=`remove_duplicates ${MANPATH}`
export LIBPATH=`remove_duplicates ${LIBPATH}`
export LIBRARYPATH=`remove_duplicates ${LIBRARYPATH}`
export SYSPATH=`remove_duplicates ${SYSPATH}`
export SHLIB_PATH=`remove_duplicates ${SHLIB_PATH}`


#--------------------------------------------------------------
# Take second environment snapshot for C-shell support.
#--------------------------------------------------------------
if ( [ "${ENV_USE_SNAPSHOTS}" == "Y" ] ) then
    printenv > ${ENV_SNAPSHOT_2}
fi

#--------------------------------------------------------------
# For User defined plot libaries, check how ACC_PLOT_PACKAGE
# is defined, if not defined then set to default "pgplot"
#--------------------------------------------------------------
if ( [ "${ACC_PLOT_PACKAGE}" == "" ] ) then
    export ACC_PLOT_PACKAGE="pgplot"
fi

#--------------------------------------------------------------
# If there are some acc dependent definitions
# (from cesr_online.bashrc), then redefine them
#--------------------------------------------------------------
if [[ -n `type -t acc_dep_define` ]]; then
	  acc_dep_define
fi

#--------------------------------------------------------------
# Set ulimits to reasonable values...
# Although not recommended, a user may override these values 
# depending on where in their .bashrc acc_vars.sh is sourced.
#--------------------------------------------------------------
ulimit -S -c 0
ulimit -S -s 10240
ulimit -S -d 25165824

#--------------------------------------------------------------
# Set Build Linker Flags
#--------------------------------------------------------------
if [ -e ${ACC_UTIL}/build_flags_config ] ; then
    . ${ACC_UTIL}/build_flags_config
fi
