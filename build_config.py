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
        'domain' : 'OFFLINE',
        #'host'   : 'acc101.lns.cornell.edu',
        'host'   : 'lnx6212.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/Gmake',
                '/trunk/src/include',
#                '/trunk/src/lattice',
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
                ##ccon_det
                ##rfintl
            ]
        }
    },
    'Linux_x86_64_intel-online' : {
        'type' : 'release',
        'platform' : 'Linux_x86_64_intel',
        'basedir' : '/gfs/cesr/online/lib',
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
                '/trunk/src/CBSM/xBSM/XbsmAnalysis',
                '/trunk/src/displays',
                '/trunk/src/logit',
                '/trunk/src/magstat',
                '/trunk/src/simcon',
                '/trunk/src/err_mon',
                '/trunk/src/fastlog',
                '/trunk/src/fbph',
                '/trunk/src/gen_log',
                '/trunk/src/newin',
                '/trunk/src/show',
                '/trunk/src/synchv',
                '/trunk/src/tao_cesr',
                '/trunk/src/vacmon',
                '/trunk/src/xscope',
		'/CESR/CESR_services/intloc'
            ]
        }
    },
    'Linux_x86_64_intel' : {
        'type' : 'packages',
        'platform' : 'Linux_x86_64_intel',
        'basedir' : '/nfs/acc/libs',
        'domain' : 'OFFLINE',
        'host' : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/packages/cfortran',
                    # A copy, perhaps, to central include directory?
                    # Makefiles need to know how to find this, so perhaps not, for now.
                '/trunk/packages/forest',
                    # gmake is all that is necessary
                ###'/trunk/packages/num_recipes/recipes_f-90',
                ###    # gmake F90="ifort" NRROOT=`pwd` F90OPTS="-Bstatic -cpp -u -check bounds -check format -warn declarations" lib
                ###    # cp -p *.mod ../modules
                ###    # cp -p librecipes_f90.a ../lib/librecipes_f90.a
                ###    # gmake NRROOT=`pwd` clean
                ###    #---------debug------------------
                ###    # gmake F90="ifort" NRROOT=`pwd` F90OPTS="-Bstatic -cpp -u -check bounds -check format -warn declarations -g" lib
                ###    # cp -p *.mod ../modules
                ###    # cp -p librecipes_f90.a ../lib/librecipes_f90_g.a
                ###    # gmake NRROOT=`pwd` clean
                '/trunk/packages/num_recipes/recipes_c-ansi',
                    # gmake -fmakefile_cesr CC="gcc -DANSI" NRROOT=`pwd` lib
                    # cp -p librecipes_c-ansi.a ../lib/librecipes_c-ansi.a
                    # gmake -fmakefile_cesr CC="gcc -DANSI" NRROOT=`pwd` clean
                    #---------debug------------------
                    # gmake -fmakefile_cesr CC="gcc -DANSI -g" NRROOT=`pwd` lib
                    # cp -p librecipes_c-ansi.a ../lib/librecipes_c-ansi_g.a
                    # gmake -fmakefile_cesr CC="gcc -DANSI" NRROOT=`pwd` clean
                '/trunk/packages/xsif',
                    # gmake
                '/trunk/packages/PGPLOT',
                    # ./makemake . linux ifort_gcc
                    # gmake
                    # gmake cpg
                    # cp -p libpgplot.a ../lib/libpgplot.a
                    # cp -p libcpgplot.a ../lib/libcpgplot.a
                    # gmake clean
                    #---------debug------------------
                    # ./makemake . linux ifort_gcc_g
                    # gmake
                    # gmake cpg
                    # cp -p libpgplot.a ../lib/libpgplot_g.a
                    # cp -p libcpgplot.a ../lib/libcpgplot_g.a
                    # gmake clean
                '/trunk/packages/gsl',
                    # ./configure --prefix `pwd`/..
                    # make
                    # make install
                '/trunk/packages/fgsl',
                    # ./configure --prefix `pwd`/.. --f90 ifort --gsl `pwd`/..
                    # make
                    # make install
                '/trunk/packages/lapack',
                    # cmake .
                    # make
                    # cp lib/* ../li
                    # Symlink library to liblapack_g.a as well.
                '/trunk/packages/lapack/LAPACK95'
                    # cd SRC
                    # make single_double_complex_dcomplex
                    # cp ../lapack95.a ../../../lib/liblapack95.a
                    # Symlink library to liblapack95_g.a as well.
                    # cp lapack95_modules/* ../../modules
                # fftw3
                    # ./configure --enable-shared --disable-dependency-tracking --enable-threads --prefix=`pwd`/..
                    # make
                    # make install
                # root
                    # This employs 32-bit Python 2.7
                    # ./configure --enable-fftw3 --with-fftw3-incdir=`pwd`/../include --with-fftw3-libdir=`pwd`/../lib --disable-python --prefix=`pwd`/.. --etcdir=`pwd`/../etc
                    # make
                    # make install
                    # User must run 'thisroot.sh' script in the directory
                    #   as part of their environment setup.
            ]
        }
    }    
}



