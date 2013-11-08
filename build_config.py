#-*-python-*-
#
# build_supervisor configuration file
#-----------------------------------------------------

offline_release_build_request = [
    #'Linux_i686_intel-offline',
    'Linux_x86_64_intel-offline' 
    ]

online_release_build_request = [
    'Linux_x86_64_intel-online'
    ]

packages_build_request = [
    'Linux_i686_intel'
    ]

dist_build_request = [
    'Linux_i686_intel'
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
    'Linux_x86_64_intel-offline' : {
        'type' : 'release',
        'platform' : 'Linux_x86_64_intel',
        'basedir' : '/nfs/acc/libs',
        'util_dir' : '/nfs/acc/libs/util',
        'domain' : 'OFFLINE',
        'host'   : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/build_system',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
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
                '/trunk/src/displays',
                '/trunk/src/logit',
                '/trunk/src/magstat',
                '/trunk/src/simcon',
		'/CESR/CESR_services/intloc',
		'/trunk/src/err_mon',
		'/trunk/src/fastlog',
		'/trunk/src/newin',
		'/trunk/src/show',
		'/trunk/src/vacmon',
		'/trunk/src/xscope',
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
		'/CESR/CESR_progs/rfnet',
		'/CESR/CESR_progs/save',
		'/CESR/CESR_progs/timing',
		'/CESR/CESR_progs/vac',
		'/CESR/CESR_libs/rf',
		'/CESR/CESR_progs/crf',
		'/CESR/CESR_progs/lrf',
		'/CESR/CESR_progs/srf',
            ]
        }
    },
    'Linux_x86_64_intel-online' : {
        'type' : 'release',
        'platform' : 'Linux_x86_64_intel',
        'basedir' : '/gfs/cesr/online/lib',
        'util_dir' : '/gfs/cesr/online/lib/util',
        'domain' : 'ONLINE',
        'host'   : 'cesr109.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/build_system',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
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
                '/trunk/src/displays',
                '/trunk/src/logit',
                '/trunk/src/magstat',
                '/trunk/src/simcon',
                '/trunk/src/err_mon',
                '/trunk/src/fastlog',
                '/trunk/src/fbph',
                '/trunk/src/gen_log',
                '/trunk/src/newin',
                '/trunk/src/save',
                '/trunk/src/show',
                '/trunk/src/synchv',
                '/trunk/src/tao_cesr',
                '/trunk/src/vacmon',
                '/trunk/src/xscope',
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
		'/CESR/CESR_progs/rfnet',
		'/CESR/CESR_progs/save',
		'/CESR/CESR_progs/timing',
		'/CESR/CESR_progs/vac',
		'/CESR/CESR_libs/rf',
		'/CESR/CESR_progs/crf',
		'/CESR/CESR_progs/lrf',
		'/CESR/CESR_progs/srf',
            ]
        }
    },
    'Linux_x86_64_intel' : {
        'type' : 'packages',
        'platform' : 'Linux_x86_64_intel',
        'basedir' : '/nfs/acc/libs',
        'util_dir' : '/nfs/acc/libs/util',
        'domain' : 'OFFLINE',
        'host' : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/packages/activemq-cpp-3.7.0',
                    # mk 
		    # mkd  
		    # Please see /trunk/packages/activemq-cpp-3.7.0/acc_build for build flags  
                '/trunk/packages/cfortran',
                    # A copy, perhaps, to central include directory?
                    # Makefiles need to know how to find this, so perhaps not, for now.
                '/trunk/packages/forest',
                    # mk
		    # mkd
                '/trunk/packages/num_recipes/recipes_c-ansi',  # Used?
                    # gmake -fmakefile_cesr CC="gcc -DANSI" NRROOT=`pwd` lib
                    # cp -p librecipes_c-ansi.a ../lib/librecipes_c-ansi.a
                    # gmake -fmakefile_cesr CC="gcc -DANSI" NRROOT=`pwd` clean
                    #  ---debug------------------
                    # gmake -fmakefile_cesr CC="gcc -DANSI -g" NRROOT=`pwd` lib
                    # cp -p librecipes_c-ansi.a ../lib/librecipes_c-ansi_g.a
                    # gmake -fmakefile_cesr CC="gcc -DANSI" NRROOT=`pwd` clean
                '/trunk/packages/xsif',
                    # mk
		    # mkd
                '/trunk/packages/PGPLOT',
                    # mk 
		    # mkd  
		    # Please see /trunk/packages/PGPLOT/acc_build for build flags  
                '/trunk/packages/gsl',
                    # mk 
		    # mkd  
		    # Please see /trunk/packages/gsl/acc_build for build flags  
                '/trunk/packages/fgsl',
                    # mk 
		    # mkd  
		    # Please see /trunk/packages/fgsl/acc_build for build flags  
                '/trunk/packages/lapack',
                    # mk 
		    # mkd  
		    # Please see /trunk/packages/lapack/acc_build for build flags  
                '/trunk/packages/lapack/LAPACK95',
		    # Gets built when lapack is built.
                    # Please see /trunk/packages/lapack/acc_build_lapack95 for build flags  
                'fftw3',
		    # Production pass
		    # ./configure --enable-shared --disable-dependency-tracking --enable-threads --prefix=`pwd`/../production --includedir=`pwd`/../include
		    # make
		    # make install

		    # Debug pass
		    # ./configure --enable-shared --disable-dependency-tracking --enable-threads --prefix=`pwd`/../debug --includedir=`pwd`/../include CFLAGS='-g -O0' FFLAGS='-g -O0'
		    # make
		    # make install
                'root',
                    # This employs 32-bit Python 2.7
		    # ./configure --enable-fftw3 --with-fftw3-incdir=`pwd`/../production/include --with-fftw3-libdir=`pwd`/../production/lib --disable-python --prefix=`pwd`/../production --etcdir=`pwd`/../production/etc --incdir=`pwd`/../include/root --enable-soversion
		    # make
		    # make install
		    # make clean
		    # ./configure --enable-fftw3 --with-fftw3-incdir=`pwd`/../debug/include --with-fftw3-libdir=`pwd`/../debug/lib --disable-python --prefix=`pwd`/../debug --etcdir=`pwd`/../debug/etc --incdir=`pwd`/../include/root --build=debug --enable-soversion
		    # make
		    # make install
                    # FOR AN ONLINE BUILD:
                    # perl -pi~ -e 's!/.fs/cesr/online!\$\{CESR_ONLINE}!' bin/thisroot.sh
                'xraylib',
                    # mk
                    # mkd
                    # Please see /trunk/packages/xraylib/acc_build for build flags

            ]
        }
    }    
}



