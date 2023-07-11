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
#offline_host = 'acc101.classe.cornell.edu'
offline_host = '$HOSTNAME'

online_base_dir = '/nfs/cesr/online/lib'
online_util_dir = online_base_dir + '/util'
online_host = 'cesr109.classe.cornell.edu'

local_base_dir = '/mnt/acc/libs'
local_util_dir = local_base_dir + '/util'
#local_host = 'lnx7179.classe.cornell.edu'
local_host = '$HOSTNAME'

makefile_dir = '/usr/bin'

release_mail_list = '$USER@cornell.edu,sbp8@cornell.edu'
local_mail_list = '$USER@cornell.edu'

#-----------------------------------------------------
#-----------------------------------------------------

release_build_list = [
                'git/bmad-doc',     #RT 60649
                '/trunk/util',
                '/trunk/build_system',
                '/trunk/include',
                '/trunk/c_utils',
                '/trunk/src/mad_tpsa',
                '/packages/forest',
                '/trunk/src/sim_utils',
                '/trunk/src/bmad',
                '/trunk/src/cpp_bmad_interface',
                '/CESR/CESR_libs/cesr_utils',
                '/CESR/CESR_libs/genplt',
                '/CESR/CESR_libs/mpmnet',
                '/CESR/CESR_libs/timing',
                '/CESR/CESR_instr/CesrBPM',
                '/CESR/CESR_instr/instr_utils',
                '/Comm/Comm_libs/cbi_net',
                '/CESR/CESR_progs/cbpmfio',
                '/CESR/CESR_instr/BeamInstSupport',
                '/CESR/CESR_instr/CBPM-TSHARC',
                '/Comm/Comm_libs/rfnet',
                '/CESR/CESR_instr/nonlin_bpm',
                '/CESR/CESR_libs/mpm_utils',
                '/CESR/CESR_libs/rf',
                '/CESR/CESR_progs/tune',
                '/CESR/CESR_progs/gen_gui',
                '/CESR/CESR_progs/diagnose',
                '/CESR/CESR_services/automail',
                '/CESR/CESR_services/averager',
                '/CESR/CESR_services/condx',
                '/CESR/CESR_services/displays',
                '/CESR/CESR_services/dt80_logger',
                '/CESR/CESR_services/err_mon',
                '/CESR/CESR_services/event_wat',
                '/CESR/CESR_services/fastlog',
                '/CESR/CESR_services/gpib_serv',
                '/CESR/CESR_services/htcmon',
                '/CESR/CESR_services/intloc',
                '/CESR/CESR_services/logit',
                '/CESR/CESR_services/onoff',
                '/CESR/CESR_services/per_mag',
                '/CESR/CESR_services/rfintl',
                '/CESR/CESR_services/sentry',
                '/CESR/CESR_services/show',
                '/CESR/CESR_services/synring',
                '/CESR/CESR_services/vacmon',
                '/CESR/CESR_services/xscope',
                '/CESR/CESR_progs/magstat',
                '/CESR/CESR_services/simcon',
                '/trunk/src/tao',
                '/trunk/src/lux',
                '/trunk/src/bmadz',
                '/trunk/src/bsim',
                '/CESR/CESR_progs/synchv',
                '/CESR/CESR_progs/cesrv',
                '/trunk/src/regression_tests',
                '/trunk/src/bsim_cesr',
                '/CESR/CESR_progs/BPM_tbt_gain',
                '/CESR/CESR_progs/cesr_programs',
                '/trunk/src/util_programs',
                '/CESR/CESR_services/CBIC',
                '/trunk/src/examples',
                '/trunk/src/analyzer',
                '/CESR/CESR_progs/xbus_book',
                '/CESR/CESR_progs/newin',
                '/CESR/CESR_progs/DB_utils',
                '/CESR/CESR_progs/chfeed',
                '/CESR/CESR_progs/gdl',
                '/CESR/CESR_progs/hard',
                '/CESR/CESR_progs/lat_utils',
                '/CESR/CESR_progs/magnet',
                '/CESR/CESR_progs/save',
                '/CESR/CESR_progs/vac',
                '/CESR/CESR_progs/crf',
                '/CESR/CESR_progs/srf',
                '/CESR/CESR_progs/univ_tune_tracker',
                '/CESR/CESR_services/console',
                '/CESR/CESR_services/winj',
                '/CESR/CESR_services/daily',
                '/CESR/CESR_services/xetec',
                '/CESR/CESR_services/webrep',
                '/CESR/CESR_services/srf232',
                '/CESR/CESR_services/bcmserv',
                '/CESR/CESR_services/moore232',
                '/CESR/CESR_services/yoko232',
                '/CESR/CESR_services/scwiggler',
                '/CESR/CESR_services/mooreenet',
                '/CESR/CESR_services/lt107_mon',
                '/CESR/CESR_services/delphi',
                '/CESR/CESR_services/runlog',
                '/CESR/CESR_services/disp_tunes',
                '/CESR/CESR_services/gen_log',
                '/CESR/CESR_services/bpm_poll',
                '/CESR/CESR_services/comet',
                '/CESR/CESR_services/powermonitor_check',
                '/CESR/CESR_progs/auto_char',
                '/CESR/CESR_progs/beam_dose',
                '/CESR/CESR_progs/beam_optimizer',
                '/CESR/CESR_progs/cbpm_mon',
                '/CESR/CESR_progs/dtp',
                '/CESR/CESR_progs/electest',
                '/CESR/CESR_progs/ethscope',
                '/CESR/CESR_progs/fbph',
                '/CESR/CESR_progs/gdl_inp',
                '/CESR/CESR_progs/grofix',
                '/CESR/CESR_progs/inj',
                '/CESR/CESR_progs/ldinit',
                '/CESR/CESR_progs/linevolt',
                '/CESR/CESR_progs/linac',
                '/CESR/CESR_services/linmon',
                '/CESR/CESR_progs/mugshot',
                '/CESR/CESR_progs/nmr_test',
                '/CESR/CESR_progs/node_set',
                '/CESR/CESR_progs/plottunes',
                '/CESR/CESR_progs/res_meas',
                '/CESR/CESR_progs/scopeget',
                '/CESR/CESR_progs/timing_test',
                '/CESR/CESR_progs/tools',
                '/CESR/CESR_progs/refresh',
                '/CESR/CESR_progs/analyze_transient',
                '/CESR/CESR_progs/sig_acq',
                '/CESR/CESR_progs/knobs',
                '/CESR/CESR_progs/xbus_load',
                '/CESR/CESR_progs/moorecon',
                '/CESR/CESR_progs/vacmap',
]

packages_build_list = [
                '/packages/recipes_f-90_LEPP',
                '/packages/activemq-cpp-3.7.0',
                '/packages/cfortran',
                '/packages/num_recipes/recipes_c-ansi',
                '/packages/xsif',
                '/packages/PGPLOT',
                '/packages/plplot',
                '/packages/gsl',
                '/packages/fgsl',
                '/packages/lapack',
                '/packages/lapack95',
                '/packages/fftw',
                '/packages/root',
                '/packages/xraylib',
                '/packages/openmpi',
                '/packages/hdf5',
                '/packages/jsonfortran',
                '/packages/libzmq',
]

#-----------------------------------------------------
#-----------------------------------------------------
repository_addresses = {
#    'ACC-CLASSE'        : 'https://accserv.classe.cornell.edu/svn',
#    'ACC-CLASSE-local'  : '/mnt/svn',
    'ACC-CLASSE'        : 'https://accserv.lepp.cornell.edu/svn',
    'ACC-CLASSE-local'  : '/mnt/svn',
    'ACC-LEPP'          : 'https://accserv.lepp.cornell.edu/svn',
    'ACC-LEPP-local'    : '/mnt/svn',
    'GitLab'            : 'https://gitlab01.classe.cornell.edu/bmad/'     #RT 60649
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
        'email_list'   : release_mail_list,
        'repositories' : {
            'ACC-CLASSE' : release_build_list
        }
    },
    'Linux_x86_64_intel-online' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : online_base_dir,
        'util_dir'     : online_util_dir,
        'domain'       : 'ONLINE',
        'host'         : online_host,
        'email_list'   : release_mail_list,
        'repositories' : {
            'ACC-CLASSE' : release_build_list
        }
    },
    'Linux_x86_64_intel-local' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'email_list'   : local_mail_list,
        'repositories' : {
            'ACC-CLASSE' : release_build_list
        }
    },
    'Linux_x86_64_gfortran-offline' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'email_list'   : release_mail_list,
        'repositories' : {
            'ACC-CLASSE' : release_build_list
        }
    },
    'Linux_x86_64_gfortran-online' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : online_base_dir,
        'util_dir'     : online_util_dir,
        'domain'       : 'ONLINE',
        'host'         : online_host,
        'email_list'   : release_mail_list,
        'repositories' : {
            'ACC-CLASSE' : release_build_list
        }
    },
    'Linux_x86_64_gfortran-local' : {
        'type'         : 'release',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'email_list'   : local_mail_list,
        'repositories' : {
            'ACC-CLASSE' : release_build_list
        }
    },
    'packages_intel'   : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'email_list'   : release_mail_list,
        'repositories' : {
            'ACC-CLASSE' : packages_build_list
        }
    },
    'packages_intel-local'   : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_intel',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'email_list'   : local_mail_list,
        'repositories' : {
            'ACC-CLASSE' : packages_build_list
        }
    },    
    'packages_gfortran' : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : offline_base_dir,
        'util_dir'     : offline_util_dir,
        'domain'       : 'OFFLINE',
        'host'         : offline_host,
        'email_list'   : release_mail_list,
        'repositories' : {
            'ACC-CLASSE' : packages_build_list
        }
    },
    'packages_gfortran-local' : {
        'type'         : 'packages',
        'platform'     : 'Linux_x86_64_gfortran',
        'basedir'      : local_base_dir,
        'util_dir'     : local_util_dir,
        'domain'       : 'LOCAL',
        'host'         : local_host,
        'email_list'   : local_mail_list,
        'repositories' : {
            'ACC-CLASSE' : packages_build_list
        }
    }    
}
