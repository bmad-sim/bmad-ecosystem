#!/usr/bin/env python

#+
# Script to test if package libraries for Bmad in the Cornell SVN repo have changed.
# If so, must port the changes to Git.
#-

import subprocess
svn_repo = {
  "https://accserv.classe.cornell.edu/svn/trunk/src/open_spacecharge"  : "50675",
  "https://accserv.classe.cornell.edu/svn/packages/PGPLOT"             : "55094",  #  Change from 50204
  "https://accserv.classe.cornell.edu/svn/packages/fftw"               : "54912",  #  Change from 54142
  "https://accserv.classe.cornell.edu/svn/packages/fgsl"               : "54914",  #  Change from 54143
  "https://accserv.classe.cornell.edu/svn/packages/gnu_utilities_src"  : "51479",
  "https://accserv.classe.cornell.edu/svn/packages/gsl"                : "54913",  #  Change from 53177
  "https://accserv.classe.cornell.edu/svn/packages/hdf5"               : "55115",  #  Change from 54917
  "https://accserv.classe.cornell.edu/svn/packages/lapack"             : "54915",  #  Change from 53213
  "https://accserv.classe.cornell.edu/svn/packages/lapack95"           : "53051",
  "https://accserv.classe.cornell.edu/svn/packages/openmpi"            : "51988",
  "https://accserv.classe.cornell.edu/svn/packages/plplot"             : "55075",  #  Change from 54925
  "https://accserv.classe.cornell.edu/svn/packages/xraylib"            : "54916",  #  Change from 53181
  "https://accserv.classe.cornell.edu/svn/trunk/util"                  : "55177",  #  Change from 55117
  "https://accserv.classe.cornell.edu/svn/trunk/build_system"          : "54615",
}

all_pass = True
for repo, rev in svn_repo.items():
  p = subprocess.Popen(f"svn info {repo} | grep 'Last Changed Rev'", stdout=subprocess.PIPE, shell=True)
  (output, err) = p.communicate()
  output = output.decode("utf-8").split(':')[1].strip()
  rep = f'"{repo}"'
  if output == rev:
    print (f'  {rep : <68} : "{output}",')
  else:
    print (f'  {rep : <68} : "{output}",  # \033[1m\033[91m Change from {rev}\033[0m')
    all_pass = False

if all_pass:
  print ('No SVN repo change.')
else:
  print ('\033[1m\033[91mSVN REPO CHANGE!!!!!!\033[0m')
