#-*-python-*-
# Accepts full build logfile name as only argument.
# i.e.
# builder /nfs/acc/libs/Linux_i686_ifort/log/cesr_2011_0907_d_1.log
#--------------------------------------------------------------
import sys
import os
from os.path import normpath, basename
import subprocess as sub
import socket

logfile = sys.argv[1]

file = open(logfile, 'r+')
inlines = file.readlines()
checkout_manifest = {}

print 'BUILDER SCRIPT RUNNING on host: ' + socket.gethostname()

#--------------------------------------------------------------
# Extract header values from log file to use as control inputs.
# Any header field found becomes a variable in the class
# 'invars' and can be accessed afterwards via the syntax
# 'invars.<var>'.
#--------------------------------------------------------------
class invars():
    pass

for line in inlines:
    if '[builder]' in line:
        break
    if 'repository' in line:
        repo = line.split()[1]
        checkout_manifest[repo] = line.split()[2:]
    else:
        var = line.split()[0]
        val = line.split()[1].strip()
        setattr(invars, var, val)



# Close file and reopen in append mode.
file.close()
file = open(logfile, 'a', 1) # unbuffered, needed?

sys.stdout = file

hostname = socket.gethostname()
p = sub.Popen('kinit -k -t /etc/cesrulib-keytab cesrulib/'+hostname,
              bufsize=1,
              shell=True,
              stdout=sub.PIPE )
while True:
    nextline = p.stdout.readline()
    if nextline == '' and p.poll() != None:
        break
    sys.stdout.write(nextline)
    sys.stdout.flush()

#print 'Shell Environment Dump:'
#for envvar in os.environ:
#    print envvar + ' = ' + os.environ[envvar]


def manifest_to_build_list( manifest ):
    """Turn list of repository check-out paths into a
       simple list of buildable directories."""
    build_list = []
    for repo in manifest:
        for dir in checkout_manifest[repo]:
            full_dir = normpath(invars.full_release_dir) + '/' + basename(normpath(dir))
            if os.path.exists(full_dir + '/Makefile'):
                print full_dir
                build_list.append(full_dir)
    return build_list


def check_out_files( manifest ):
    """Check out all files described in the manifest
       associated with this build."""
    print 'Checking out files according to manifest...'
    svnlines = []
    for repo in manifest:
        svnlines = []
        for dir in checkout_manifest[repo]:
            p = sub.Popen('kinit -k -t /etc/cesrulib-keytab cesrulib/'+hostname+'; svn co ' + repo + dir,
                          bufsize=1,
                          shell=True,
                          stdout=sub.PIPE,
                          stderr=sub.STDOUT )
            while True:
                nextline = p.stdout.readline()
                if nextline == '' and p.poll() != None:
                    break
                sys.stdout.write(nextline)
                sys.stdout.flush()


def link_to_packages( packages_name ):
    """Create a symbolic link in the release directory
       called 'packages' to the packages area named in
       the argument."""
    packages_dir = invars.libs_basedir+'/'+invars.platform+'/'+packages_name
    os.symlink( packages_dir, invars.full_release_dir+'/packages' )


#def determine_build_order( ):
#"""Examine all source code to build and come up with
#   optimum order of directories to visit."""


def build_directory( dir ):
    print '\n\n\n-------- Building: ' + dir
    os.chdir( dir )
    build_command = 'ACCLIB='+invars.build_name+'; ACC_FORCE_32_BIT=N; source ' +invars.util_dir+'/acc_vars.sh; ifort -v; printenv | grep ACC; gmake PRECISION="_DBL" DO_EXTRA_MAKES=Y USE_PGPLOT=Y'
    p = sub.Popen(build_command,
                  bufsize=1,
                  shell=True,
                  stdout=sub.PIPE,
                  stderr=sub.STDOUT )
    while True:
        nextline = p.stdout.readline()
        if nextline == '' and p.poll() != None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()


#---------------------------
        
link_to_packages( invars.packages_name )


blist = manifest_to_build_list( checkout_manifest )
for dir in blist:
    build_directory( dir )

    
