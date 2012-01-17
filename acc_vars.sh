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

#HAC_ARCHIVE_BASE_DIR='/gfs/cesr/online/lib'
#HAC_RELEASE_MGMT_DIR='/gfs/cesr/online/lib/util'
#HAC_OPT_SOFTWARE_DIR='/nfs/opt'

HAC_ARCHIVE_BASE_DIR='/nfs/acc/libs'
HAC_RELEASE_MGMT_DIR='/nfs/acc/libs/util'
HAC_OPT_SOFTWARE_DIR='/nfs/opt'

#ONLINE_ARCHIVE_BASE_DIR='/nfs/cesr/online/lib'
#ONLINE_RELEASE_MGMT_DIR='/nfs/cesr/online/lib/util'
#ONLINE_OPT_SOFTWARE_DIR='/nfs/opt'

ONLINE_ARCHIVE_BASE_DIR='/nfs/acc/libs'
ONLINE_RELEASE_MGMT_DIR='/nfs/acc/libs/util'
ONLINE_OPT_SOFTWARE_DIR='/nfs/opt'

OFFLINE_ARCHIVE_BASE_DIR='/nfs/acc/libs'
OFFLINE_RELEASE_MGMT_DIR='/nfs/acc/libs/util'
OFFLINE_OPT_SOFTWARE_DIR='/nfs/opt'


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




#--------------------------------------------------------------
# Support functions
#--------------------------------------------------------------
# Is argument $1 missing from argument $2 (or PATH)?
function no_path() {
        eval "case :\$${2-PATH}: in *:$1:*) return 1;; *) return 0;; esac"
}

# If argument path ($1) exists and is not in path, append it.
function add_path () {
  [ -d ${1:-.} ] && no_path $* && eval ${2:-PATH}="\$${2:-PATH}:$1"
}

# If argument ($1) is in path, remove it.
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
# Set CESR_ONLINE to the path used for 'CESR offline'
# hosts if it hasn't already been set upon entering this
# script.
#--------------------------------------------------------------
CESR_ONLINE=${CESR_ONLINE:-/nfs/cesr/online}


#--------------------------------------------------------------
# 'CESR ONLINE' hosts are interactive shell capable machines in
# that follow the convention of having 'cesr' in their hostname.
# If the machine on which this script is sourced is defined to
# be a 'CESR_ONLINE' host, the follow value will be 
# shell-TRUE(0) and (1) otherwise.
#
#  TODO: Double check that no laboratory hostnames violate
#        this convention.  'cesrweb'?
#--------------------------------------------------------------
if ( [ `hostname | grep cesr` ] ) then
    IS_CESR_ONLINE_HOST='true'
fi


#--- CESR ONLINE Hosts
if ( [ "${IS_CESR_ONLINE_HOST}" == "true" ] ) then

    #-- CESR HAC (High-availability Cluster) Online Host
    if ( [ "${CESR_ONLINE}" == "/gfs/cesr/online" ] )then
	SETUP_SCRIPTS_DIR=${HAC_RELEASE_MGMT_DIR}
	RELEASE_ARCHIVE_BASE_DIR=${HAC_ARCHIVE_BASE_DIR}
	OPT_SOFTWARE_DIR=${HAC_OPT_SOFTWARE_DIR}
    else
    #-- CESR Online Host
	SETUP_SCRIPTS_DIR=${ONLINE_RELEASE_MGMT_DIR}
	RELEASE_ARCHIVE_BASE_DIR=${ONLINE_ARCHIVE_BASE_DIR}
	OPT_SOFTWARE_DIR=${ONLINE_OPT_SOFTWARE_DIR}
    fi

#-- CESR OFFLINE Host
else
    SETUP_SCRIPTS_DIR=${OFFLINE_RELEASE_MGMT_DIR}
    RELEASE_ARCHIVE_BASE_DIR=${OFFLINE_ARCHIVE_BASE_DIR}
    OPT_SOFTWARE_DIR=${OFFLINE_OPT_SOFTWARE_DIR}

fi 


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

export ACC_GMAKE=${ACC_RELEASE_DIR}/Gmake
export CESR_GMAKE=${ACC_GMAKE}  # For backwards compatibility.



#--------------------------------------------------------------
# Non-ACC naming convention variables. 
# Named such for historical compatibility.
# FIXME:  CESR_ONLINE  should be moved to a different script.
#--------------------------------------------------------------
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
# architecture.
#
# TODO: These vendor-supplied compiler setup scripts will
#       need to live in both CESR OFFLINE and CESR ONLINE 
#       directories.
#--------------------------------------------------------------
case ${ACC_OS_ARCH} in

    "Linux_i686" )
        # 32-bit Intel Fortran (ifort)
	source ${OPT_SOFTWARE_DIR}/ifc/bin/ifortvars.sh
	;;

    "Linux_x86_64" )
	# Add ifort-9-specific path information to user's environment
	# to allow running 32-bit fortran programs on a 64-bit machine.
	# Then source the final 64-bit path information on top.
	source ${OPT_SOFTWARE_DIR}/ifc/bin/ifortvars.sh
        # 64-bit Intel Fortran (ifort) v12.1.0.233
	source ${OPT_SOFTWARE_DIR}/intel/composerxe/bin/compilervars.sh intel64
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
# Ensure that there is only one reference to the release's
# solib directory by searching the LD_LIBRARY_PATH variable
# value for paths that are prefixed with the
# RELEASE_ARCHIVE_BASE_DIR and removing them.  The last
# step is appending the presently valid solib path for
# the active release.
#--------------------------------------------------------------
LD_LIBRARY_PATH_TEMP=${LD_LIBRARY_PATH}
LD_LIBRARY_PATH=""
DEFAULT_IFS=${IFS} # bash 'internal field separator' reserved variable
IFS=":"
for part in ${LD_LIBRARY_PATH_TEMP}; do
  case ${part} in

      "${RELEASE_ARCHIVE_BASE_DIR}"*)
	  #echo "ACC path found. Will remove from LD_LIBRARY_PATH."
	  continue
	  ;;

      *)
	  if ( [ "${LD_LIBRARY_PATH}" == "" ] ) then
	      LD_LIBRARY_PATH=${part}
	  else
	      LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${part}
	  fi  

  esac
done
IFS=${DEFAULT_IFS}

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ACC_RELEASE_DIR}/solib



#--------------------------------------------------------------
# Variables to locate and use CERNlib, ROOT.
# To be integrated into packages system for CESR management?
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


alias current='ACC_RELEASE_REQUEST=current; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'
alias devel='ACC_RELEASE_REQUEST=devel; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh'

alias current64='ACC_FORCE_32_BIT=N; current'
alias devel64='ACC_FORCE_32_BIT=N; devel'

# These are functions so that they may take arguments from the user.

function setrel() {
    ACC_RELEASE_REQUEST=${1}; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh
}

function setrel64() {
    ACC_FORCE_32_BIT=N; ACC_RELEASE_REQUEST=${1}; source ${SETUP_SCRIPTS_DIR}/acc_vars.sh
}



#--------------------------------------------------------------
# Take second environment snapshot for C-shell support.
#--------------------------------------------------------------
if ( [ "${ENV_USE_SNAPSHOTS}" == "Y" ] ) then
    printenv > ${ENV_SNAPSHOT_2}
fi
