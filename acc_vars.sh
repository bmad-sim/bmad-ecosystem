#--------------------------------------------------------------
# acc_vars_sh
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
# For C-shell variant use, all of the above applies, but the
# user must source the script named 
#   acc_vars.csh  instead.
#--------------------------------------------------------------

#--------------------------------------------------------------
# Files to hold a 'printenv' environment snapshot before this
# script defines all its variables, and again after for use
# in the wrapper script that C-style shells only will source.
# This is so that all the syntax remains in one language,
# (bash) and need not be kept up-to-date between a bourne-shell
# setup script and C-shell one.
#
# Take first environment snapshot for C-shell support.
#--------------------------------------------------------------

# Capture value of ACC_BIN to allow removal from path for cleanliness.
OLD_ACC_BIN=${ACC_BIN}

# Temporary for cleanliness.
source /home/cesrulib/bin/clean_env.sh


# Change this to /nfs/acc/libs/bin once login scripts get finalized
USER_BIN_DIR=/home/cesrulib/bin


#--------------------------------------------------------------
# Support functions
#--------------------------------------------------------------
# is $1 missing from $2 (or PATH) ?
function no_path() {
        eval "case :\$${2-PATH}: in *:$1:*) return 1;; *) return 0;; esac"
}

# if $1 exists and is not in path, append it
function add_path () {
  [ -d ${1:-.} ] && no_path $* && eval ${2:-PATH}="\$${2:-PATH}:$1"
}

# if $1 is in path, remove it
function del_path () {
  no_path $* || eval ${2:-PATH}=`eval echo :'$'${2:-PATH}: |
    /bin/sed -e "s;:$1:;:;g" -e "s;^:;;" -e "s;:\$;;"`
}


if ( [ "${ENV_USE_SNAPSHOTS}" == "Y" ] ) then
    ENV_SNAPSHOT_1="${HOME}/.ACC_env_snapshot_1.tmp"
    ENV_SNAPSHOT_2="${HOME}/.ACC_env_snapshot_2.tmp"
    printenv > ${ENV_SNAPSHOT_1}
fi



#--------------------------------------------------------------
# Set the location of the directory containing all the 
# necessary environment setup scripts.
#--------------------------------------------------------------
SETUP_SCRIPTS_DIR=/nfs/acc/libs/util



#--------------------------------------------------------------
# Central area where all releases are kept.
#--------------------------------------------------------------
RELEASE_ARCHIVE_BASE_DIR=/nfs/acc/libs



#--------------------------------------------------------------
# Variables that influence the build process
#--------------------------------------------------------------
ACC_OS="`uname`"
export ACC_ARCH="`uname -m`"
export ACC_OS_ARCH="${ACC_OS}_${ACC_ARCH}"
case ${ACC_OS_ARCH} in

    "Linux_i686" )
	export ACC_FC="intel";;

    "Linux_x86_64" )
	export ACC_FC="intel";;

    "Darwin_i386" )
	export ACC_FC="gfortran";;

    * )
	export ACC_FC="gfortran";;
esac



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
if ( [ "${ACC_FORCE_32_BIT}" == "" ] ) then
    ACC_FORCE_32_BIT=Y
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
export ACC_TRUE_RELEASE=`readlink ${ACC_RELEASE_DIR}`
if ( [ "${ACC_TRUE_RELEASE}" == "" ] ) then
    ACC_TRUE_RELEASE=${ACC_RELEASE}
fi
export ACC_UTIL=${SETUP_SCRIPTS_DIR}
export ACC_BIN=${ACC_RELEASE_DIR}/bin
export ACC_EXE=${ACC_BIN} # For backwards compatibility
export ACC_PKG=${ACC_RELEASE_DIR}/packages
export ACC_REPO=https://accserv.lepp.cornell.edu/svn/
export ACCR=https://accserv.lepp.cornell.edu/svn/

export ACC_GMAKE=/home/cesrulib/bin/Gmake
export CESR_GMAKE=${ACC_GMAKE}  # For backwards compatibility.



#--------------------------------------------------------------
# Non-ACC naming convention variables. 
# Named such for historical compatibility.
# FIXME:  CESR_ONLINE  should be moved to a different script.
#--------------------------------------------------------------
export CESR_ONLINE=/nfs/cesr/online
export BMAD_LAT=${ACC_RELEASE_DIR}/config/bmad
export CESR_LAT=${ACC_RELEASE_DIR}/config/cesr
export CTA_LAT=${ACC_RELEASE_DIR}/config/cta
export ERL_LAT=${ACC_RELEASE_DIR}/config/erl
export ILC_LAT=${ACC_RELEASE_DIR}/config/ilc



#--------------------------------------------------------------
# Variables to support other software packages
#--------------------------------------------------------------
export PGPLOT_DIR=${ACC_PKG}/PGPLOT
export PGPLOT_FONTS=${ACC_PKG}/PGPLOT



#--------------------------------------------------------------
# Set up the fortran compiler appropriate to the machine
# architecture
#--------------------------------------------------------------
case ${ACC_OS_ARCH} in

    "Linux_i686" )
        # 32-bit Intel Fortran (ifort)
        source /nfs/opt/ifc/bin/ifortvars.sh;;

    "Linux_x86_64" )
        # 64-bit Intel Fortran (ifort)
        #source /nfs/opt/intel/composerxe-2011.4.191/bin/compilervars.sh intel64;;
	source /nfs/opt/intel/composerxe/bin/compilervars.sh intel64;; # ifort v12.1.0.233

esac



#--------------------------------------------------------------
# Add the directory where utilities for the user are to be
# provided.  getf, listf, etc...
#--------------------------------------------------------------
add_path ${USER_BIN_DIR}
del_path ${OLD_ACC_BIN}
add_path ${ACC_BIN}


#--------------------------------------------------------------
# Variables to locate and use CERNlib
#--------------------------------------------------------------
export CERN_ROOT=/nfs/cern/pro
export CERNLIB=${CERN_ROOT}/lib
export CERNSRC=${CERN_ROOT}/src


#--------------------------------------------------------------
# Useful aliases
#--------------------------------------------------------------
# Information about build system setup for the current shell session.
alias accinfo=' echo -e "ACC release      : \"${ACC_RELEASE}\" [${ACC_TRUE_RELEASE}]
Architecture     : ${ACC_ARCH}
Fortran compiler : ${ACC_FC}" '
alias ACCINFO='accinfo'
alias acc_info='accinfo'
alias ACC_INFO='accinfo'



#--------------------------------------------------------------
# Take second environment snapshot for C-shell support.
#--------------------------------------------------------------
if ( [ "${ENV_USE_SNAPSHOTS}" == "Y" ] ) then
    printenv > ${ENV_SNAPSHOT_2}
fi
