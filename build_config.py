#-*-python-*-
#
# build_supervisor configuration file
#-----------------------------------------------------

intel_offline_release_build_request = [
    'Linux_x86_64_intel-offline'
    ]

intel_online_release_build_request = [
    'Linux_x86_64_intel-online'
    ]

intel_packages_build_request = [
    'packages_intel'
    ]

intel_dist_build_request = [
    'Linux_i686_intel'
    ]

intel_local_release_build_request = [
    'Linux_x86_64_intel-local'
    ]

intel_local_packages_build_request = [
    'packages_intel-local'
    ]

gfortran_offline_release_build_request = [
    'Linux_x86_64_gfortran-offline' 
    ]

gfortran_online_release_build_request = [
    'Linux_x86_64_gfortran-online'
    ]

gfortran_packages_build_request = [
    'packages_gfortran'
    ]

gfortran_dist_build_request = [
    'Linux_i686_gfortran'
    ]

gfortran_local_release_build_request = [
    'Linux_x86_64_gfortran-local'
    ]

gfortran_local_packages_build_request = [
    'packages_gfortran-local'
    ]


#-----------------------------------------------------
# Collect all build requests by type into a master
# dictionary.
#-----------------------------------------------------
build_requests = {}
build_requests['release_intel'] = intel_offline_release_build_request
build_requests['online-release_intel'] = intel_online_release_build_request
build_requests['packages_intel'] = intel_packages_build_request
build_requests['dist_intel'] = intel_dist_build_request
build_requests['local-release_intel'] = intel_local_release_build_request
build_requests['local-packages_intel'] = intel_local_packages_build_request

build_requests['release_gfortran'] = gfortran_offline_release_build_request
build_requests['online-release_gfortran'] = gfortran_online_release_build_request
build_requests['packages_gfortran'] = gfortran_packages_build_request
build_requests['dist_gfortran'] = gfortran_dist_build_request
build_requests['local-release_gfortran'] = gfortran_local_release_build_request
build_requests['local-packages_gfortran'] = gfortran_local_packages_build_request

#-----------------------------------------------------
#-----------------------------------------------------
offline_base_dir = '/nfs/acc/libs'
offline_util_dir = offline_base_dir + '/util'
offline_host = 'acc101.lns.cornell.edu'

online_base_dir = '/gfs/cesr/online/lib'
online_util_dir = online_base_dir + '/util'
online_host = 'cesr109.lns.cornell.edu'

local_base_dir = '/mnt/acc/libs'
local_util_dir = local_base_dir + '/util'
local_host = 'lnx7179.lns.cornell.edu'

makefile_dir = '/home/cesrulib/bin/Gmake'


#-----------------------------------------------------
#-----------------------------------------------------

local_build_list = [
                '/trunk/util',
                '/trunk/build_system',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
                '/CESR/CESR_libs/timing',
                '/trunk/src/bmad',
                '/trunk/src/cesr_utils',
                '/CESR/CESR_instr/instr_utils',
                '/trunk/src/cbi_net',
                '/trunk/src/cbpmfio',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/Comm/Comm_libs/rfnet',
                '/trunk/src/CesrBPM',                
                '/trunk/src/mpm_utils',
                '/trunk/src/nonlin_bpm',
                '/trunk/src/tao',
                '/trunk/src/tao_cesr',
                '/trunk/src/bmadz',
                '/trunk/src/bsim',
                '/trunk/src/cesrv',
                '/trunk/src/regression_tests',
                '/trunk/src/bsim_cesr',
                '/trunk/src/genplt',
                '/trunk/src/cesr_programs',
                '/trunk/src/util_programs',
                '/trunk/src/CBIC',
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples',
                '/CESR/CESR_progs/xbus_book',
                '/trunk/src/CBSM/xBSM/XbsmAnalysis',
                '/CESR/CESR_services/displays',
                '/CESR/CESR_services/logit',
                '/trunk/src/magstat',
                '/trunk/src/simcon',
                '/CESR/CESR_services/intloc',
                '/CESR/CESR_services/err_mon',
                '/CESR/CESR_services/fastlog',
                '/trunk/src/newin',
                '/CESR/CESR_services/show',
                '/CESR/CESR_services/vacmon',
                '/CESR/CESR_services/xscope',
                '/CESR/CESR_services/condx',
                '/CESR/CESR_services/dt80_logger',
                '/CESR/CESR_services/event_wat',
                '/CESR/CESR_services/htcmon',
                '/CESR/CESR_progs/DB_utils',
                '/CESR/CESR_progs/chfeed',
                '/CESR/CESR_progs/diagnose',
                '/CESR/CESR_progs/gdl',
                '/CESR/CESR_progs/hard',
                '/CESR/CESR_progs/lat_utils',
                '/CESR/CESR_progs/magnet',
                '/Comm/Comm_libs/rfnet',
                '/CESR/CESR_progs/save',
                '/CESR/CESR_progs/vac',
                '/CESR/CESR_libs/rf',
                '/CESR/CESR_progs/crf',
                '/CESR/CESR_progs/lrf',
                '/CESR/CESR_progs/srf',
]

packages_build_list = [
                '/trunk/packages/activemq-cpp-3.7.0',
                '/trunk/packages/cfortran',
                '/trunk/packages/forest',
                '/trunk/packages/num_recipes/recipes_c-ansi',
                '/trunk/packages/xsif',
                '/trunk/packages/PGPLOT',
                '/trunk/packages/plplot',
                '/trunk/packages/gsl',
                '/trunk/packages/fgsl',
                '/trunk/packages/lapack',
                '/trunk/packages/fftw',
                '/trunk/packages/root',
                '/trunk/packages/xraylib',
]

#-----------------------------------------------------
#-----------------------------------------------------
repository_addresses = {
    'ACC-LEPP'        : 'https://accserv.lepp.cornell.edu/svn',
    'ACC-LEPP-local'  : '/mnt/svn',
    'UAP-Sourceforge' : 'https://accelerator-ml.svn.sourceforge.net/svnroot/accelerator-ml/uap'
    }


#-----------------------------------------------------
#-----------------------------------------------------
build_specs = {
    'Linux_x86_64_intel-offline' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'repositories' : {
            'ACC-LEPP' : local_build_list
        }
    },
    'Linux_x86_64_intel-online' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : online_base_dir,
        'util_dir'     : online_util_dir,
        'domain'       : 'ONLINE',
        'host'         : online_host,
        'repositories' : {
            'ACC-LEPP' : local_build_list
        }
    },
    'Linux_x86_64_intel-local' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'repositories' : {
            'ACC-LEPP' : local_build_list
        }
    },
    'Linux_x86_64_gfortran-offline' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'repositories' : {
            'ACC-LEPP' : local_build_list
        }
    },
    'Linux_x86_64_gfortran-online' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : online_base_dir,
        'util_dir'     : online_util_dir,
        'domain'       : 'ONLINE',
        'host'         : online_host,
        'repositories' : {
            'ACC-LEPP' : local_build_list
        }
    },
    'Linux_x86_64_gfortran-local' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'repositories' : {
            'ACC-LEPP' : local_build_list
        }
    },
    'packages_intel'   : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'repositories' : {
            'ACC-LEPP' : packages_build_list
        }
    },
    'packages_intel-local'   : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'repositories' : {
            'ACC-LEPP' : packages_build_list
        }
    },    
    'packages_gfortran' : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'repositories' : {
            'ACC-LEPP' : packages_build_list
        }
    },
    'packages_gfortran-local' : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'repositories' : {
            'ACC-LEPP' : packages_build_list
        }
    }    
}
