#
#-*-python-*-
#-----------------------------------------------------

release_build_request = [
    'Linux_i686_intel',
    'Linux_x86_64_intel'
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
build_requests['RELEASE'] = release_build_request
build_requests['PACKAGES'] = packages_build_request
build_requests['DIST'] = dist_build_request


#-----------------------------------------------------
#-----------------------------------------------------
libs_basedir = '/nfs/acc/libs'
util_dir = '/nfs/acc/libs/util'
makefile_dir = '/home/cesrulib/bin/Gmake-testing'


#-----------------------------------------------------
#-----------------------------------------------------
repository_addresses = {
    'ACC-LEPP' : 'https://accserv.lepp.cornell.edu/svn',
    'ACC-LEPP-local' : '/mnt/svn',
    'UAP-Sourceforge' : 'https://accelerator-ml.svn.sourceforge.net/svnroot/accelerator-ml/uap'
    }


#-----------------------------------------------------
#-----------------------------------------------------
release_build_specs = {
    'Linux_i686_intel' : {
        'host' : 'lnx209.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/Gmake',
                '/trunk/src/include',
                '/trunk/src/lattice',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/CesrBPM',
                '/trunk/src/mpmnet',
                '/trunk/src/cbi_net',
                '/trunk/src/bmad',
                '/trunk/src/cesr_utils',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/CBIC',
                '/trunk/src/mpm_utils',
                '/trunk/src/nonlin_bpm',
                '/trunk/src/tao',
                '/trunk/src/tao_cesr',
                '/trunk/src/cesr_programs',
                '/trunk/src/bmadz',
                '/trunk/src/cesrv',
                '/trunk/src/bsim',
                '/trunk/src/bsim_cesr',
                '/trunk/src/util_programs',
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples'
            ]
        }
    },
    'Linux_x86_64_intel' : {
        'host' : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/Gmake',
                '/trunk/src/include',
                '/trunk/src/lattice',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/CesrBPM',
                '/trunk/src/mpmnet',
                '/trunk/src/cbi_net',
                '/trunk/src/bmad',
                '/trunk/src/cesr_utils',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/CBIC',
                '/trunk/src/mpm_utils',
                '/trunk/src/nonlin_bpm',
                '/trunk/src/tao',
                '/trunk/src/tao_cesr',
                '/trunk/src/cesr_programs',
                '/trunk/src/bmadz',
                '/trunk/src/cesrv',
                '/trunk/src/bsim',
                '/trunk/src/bsim_cesr',
                '/trunk/src/util_programs',
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples'
            ]
        }
    },
    'OSF_alpha_hp' : {
        'host' : 'cesr66.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/Gmake',
                '/trunk/src/lattice',
                '/trunk/src/c_utils',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/sim_utils',
                '/trunk/src/mpmnet',
                '/trunk/src/cbi_net',
                '/trunk/src/BeamInstSupport',
                '/trunk/src/CBPM-TSHARC',
                '/trunk/src/bmad',
                '/trunk/src/cesr_utils',
                '/trunk/src/mpm_utils',
                '/trunk/src/nonlin_bpm',
                '/trunk/src/tao',
                '/trunk/src/tao_cesr',
                '/trunk/src/cesr_programs',
                '/trunk/src/CesrBPM',
                '/trunk/src/bmadz',
                '/trunk/src/bsim',
                '/trunk/src/bsim_cesr',
                '/trunk/src/util_programs',
                '/trunk/src/BPM_tbt_gain',
                '/trunk/src/examples'
            ]
        }
    }
}



packages_build_specs = {
    'Linux_i686_intel' : {
        'host' : 'lnx209.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/packages/cfortran',
                '/trunk/packages/forest',
                '/trunk/packages/num_recipes/recipes_f-90',
                '/trunk/packages/num_recipes/recipes_c-ansi',
                '/trunk/packages/xsif',
                '/trunk/packages/PGPLOT'
            ]
        }
    },
    'Linux_x86_64_intel' : {
        'host' : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/packages/cfortran',
                '/trunk/packages/forest',
                '/trunk/packages/num_recipes/recipes_f-90',
                '/trunk/packages/num_recipes/recipes_c-ansi',
                '/trunk/packages/xsif',
                '/trunk/packages/PGPLOT'
            ]
        }
    }    
}



dist_build_specs = {
    'Linux_i686_intel' : {
        'host' : 'lnx209.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/Gmake',
                '/trunk/src/include',
                '/trunk/src/bmad',
                '/trunk/src/examples',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/cesr_utils',
                '/trunk/packages/forest',
                '/trunk/packages/PGPLOT',
                '/trunk/packages/xsif',
                '/trunk/src/lattice',
                '/trunk/src/examples',
                '/trunk/src/sim_utils',
                '/trunk/src/bmadz',
                '/trunk/src/tao',
                '/trunk/src/bsim'
            ]
        }
    },
    'Linux_x86_64_intel' : {
        'host' : 'acc101.lns.cornell.edu',
        'repositories' : {
            'ACC-LEPP' : [
                '/trunk/util',
                '/trunk/Gmake',
                '/trunk/src/include',
                '/trunk/src/bmad',
                '/trunk/src/examples',
                '/trunk/src/recipes_f-90_LEPP',
                '/trunk/src/cesr_utils',
                '/trunk/packages/forest',
                '/trunk/packages/PGPLOT',
                '/trunk/packages/xsif',
                '/trunk/src/lattice',
                '/trunk/src/examples',
                '/trunk/src/sim_utils',
                '/trunk/src/bmadz',
                '/trunk/src/tao',
                '/trunk/src/bsim'
            ]
        }
    }
}


#-----------------------------------------------------
# Collect all platform build specifications into one
# master dictionary.
#-----------------------------------------------------
build_specs = {}
build_specs['RELEASE'] = release_build_specs
build_specs['PACKAGES'] = packages_build_specs
build_specs['DIST'] = dist_build_specs

