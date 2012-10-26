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
                '/trunk/Gmake',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
                '/trunk/src/cbi_net',
                '/trunk/src/cbpmfio',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/CBIC',
                '/trunk/src/bmad',
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
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples',
                '/trunk/src/genplt',
                '/trunk/src/CBSM/xBSM/XbsmAnalysis',
                '/trunk/src/displays',
                '/trunk/src/logit',
                '/trunk/src/magstat',
                '/trunk/src/simcon',
		'/CESR/CESR_services/intloc'
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
                '/trunk/Gmake',
                '/trunk/src/include',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
                '/trunk/src/cbi_net',
                '/trunk/src/cbpmfio',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/CBIC',
                '/trunk/src/bmad',
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
		'/CESR/CESR_services/intloc',
		'/CESR/CESR_progs/vac',
		'/CESR/CESR_progs/linac'
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
                    # ./makemake . linux ifort_gcc
                    # gmake
                    # gmake cpg
                    # cp -p libpgplot.a ../production/lib/libpgplot.a
                    # cp -p libcpgplot.a ../production/lib/libcpgplot.a
                    # gmake clean
                    #---------debug------------------
                    # ./makemake . linux ifort_gcc_g
                    # gmake
                    # gmake cpg
                    # cp -p libpgplot.a ../debug/lib/libpgplot.a
                    # cp -p libcpgplot.a ../debug/lib/libcpgplot.a
                    # gmake clean
                '/trunk/packages/gsl',
                    # ./configure --prefix `pwd`/../production
                    # make
                    # make install
                '/trunk/packages/fgsl',
                    # ./configure --prefix `pwd`/../production --f90 ifort --gsl `pwd`/../production
                    # make
                    # make install
                '/trunk/packages/lapack',
                    # cmake .
                    # make
                    # cp lib/* ../production/lib
		    # cp lib/* ../debug/lib
                '/trunk/packages/lapack/LAPACK95',
                    # cd SRC
                    # make single_double_complex_dcomplex
		    # cd ..
                    # cp lapack95.a ../../../production/lib
		    # cp lapack95.a ../../../debug/lib
                    # cp lapack95_modules/* ../../production/modules
                    # cp lapack95_modules/* ../../debug/modules
                'fftw3'
                    # ./configure --enable-shared --disable-dependency-tracking --enable-threads --prefix=`pwd`/../production
                    # make
                    # make install
                'root',
                    # This employs 32-bit Python 2.7
                    # ./configure --enable-fftw3 --with-fftw3-incdir=`pwd`/../production/include --with-fftw3-libdir=`pwd`/../production/lib --disable-python --prefix=`pwd`/../production --etcdir=`pwd`/../production/etc
                    # make
                    # make install
                    # FOR AN ONLINE BUILD:
                    # perl -pi~ -e 's!/.fs/cesr/online!\$\{CESR_ONLINE}!' bin/thisroot.sh

            ]
        }
    }    
}



