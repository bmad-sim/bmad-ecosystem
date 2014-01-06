#-*-python-*-
#
# build_supervisor configuration file
#-----------------------------------------------------

import os

offline_release_build_request = [
    'Linux_x86_64_' + os.environ["ACC_FC"] + '-offline' 
    ]

online_release_build_request = [
    'Linux_x86_64_' + os.environ["ACC_FC"] + '-online'
    ]

packages_build_request = [
    'packages_' + os.environ["ACC_FC"]
    ]

dist_build_request = [
    'Linux_i686_' + os.environ["ACC_FC"]
    ]


#-----------------------------------------------------
# Collect all build requests by type into a master
# dictionary.
#-----------------------------------------------------
build_requests = {}
build_requests['release'] = offline_release_build_request
build_requests['online-release'] = online_release_build_request
build_requests['packages'] = packages_build_request
build_requests['dist'] = dist_build_request


#-----------------------------------------------------
#-----------------------------------------------------
util_dir = '/nfs/acc/libs/util'
makefile_dir = '/home/cesrulib/bin/Gmake'


#-----------------------------------------------------
#-----------------------------------------------------
repository_addresses = {
    'ACC-LEPP' : 'https://accserv.lepp.cornell.edu/svn',
    'ACC-LEPP-local' : '/mnt/svn',
    'UAP-Sourceforge' : 'https://accelerator-ml.svn.sourceforge.net/svnroot/accelerator-ml/uap'
    }


#-----------------------------------------------------
#-----------------------------------------------------
build_specs = {
    'Linux_x86_64_' + os.environ["ACC_FC"] + '-offline' : {
        'type'     : 'release',
        'platform' : 'Linux_x86_64_' + os.environ["ACC_FC"],
        'basedir'  : '/nfs/acc/libs',
        'util_dir' : '/nfs/acc/libs/util',
        'domain'   : 'OFFLINE',
        'host'     : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/build_system',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
                '/CESR/CESR_libs/timing',
                '/CESR/CESR_instr/instr_utils',
                '/trunk/src/cbi_net',
                '/trunk/src/cbpmfio',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/bmad',
                '/Comm/Comm_libs/rfnet',
                '/trunk/src/cesr_utils',
                '/trunk/src/CesrBPM',                
                '/trunk/src/mpm_utils',
                '/trunk/src/nonlin_bpm',
                '/trunk/src/tao',
                '/trunk/src/tao_cesr',
                '/trunk/src/bmadz',
                '/trunk/src/cesrv',
                '/trunk/src/bsim',
                '/trunk/src/bsim_cesr',
                '/trunk/src/cesr_programs',
                '/trunk/src/util_programs',
                '/trunk/src/CBIC',
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples',
                '/trunk/src/genplt',
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
        }
    },
    'Linux_x86_64_' + os.environ["ACC_FC"] + '-online' : {
        'type'     : 'release',
        'platform' : 'Linux_x86_64_' + os.environ["ACC_FC"],
        'basedir'  : '/gfs/cesr/online/lib',
        'util_dir' : '/gfs/cesr/online/lib/util',
        'domain'   : 'ONLINE',
        'host'     : 'cesr109.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/build_system',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
                '/CESR/CESR_libs/timing',
                '/CESR/CESR_instr/instr_utils',
                '/trunk/src/cbi_net',
                '/trunk/src/cbpmfio',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/bmad',
                '/Comm/Comm_libs/rfnet',
                '/trunk/src/cesr_utils',
                '/trunk/src/CesrBPM',                
                '/trunk/src/mpm_utils',
                '/trunk/src/nonlin_bpm',
                '/trunk/src/tao',
                '/trunk/src/tao_cesr',
                '/trunk/src/bmadz',
                '/trunk/src/cesrv',
                '/trunk/src/bsim',
                '/trunk/src/bsim_cesr',
                '/trunk/src/cesr_programs',
                '/trunk/src/util_programs',
                '/trunk/src/CBIC',
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples',
                '/trunk/src/genplt',
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
        }
    },
    'packages_' + os.environ["ACC_FC"] : {
        'type'     : 'packages',
        'platform' : 'packages_' + os.environ["ACC_FC"],
        'basedir'  : '/nfs/acc/libs',
        'util_dir' : '/nfs/acc/libs/util',
        'domain'   : 'OFFLINE',
        'host'     : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/packages/activemq-cpp-3.7.0',
                '/trunk/packages/cfortran',
                '/trunk/packages/forest',
                '/trunk/packages/num_recipes/recipes_c-ansi',
                '/trunk/packages/xsif',
                '/trunk/packages/PGPLOT',
                '/trunk/packages/gsl',
                '/trunk/packages/fgsl',
                '/trunk/packages/lapack',
                '/trunk/packages/fftw',
                '/trunk/packages/root',
                '/trunk/packages/xraylib',
            ]
        }
    }    
}
