#!/usr/bin/env python

#+
# Script to test if package libraries for Bmad in the Cornell SVN repo have changed.
# If so, must port the changes to Git.
#-

import subprocess
svn_repo = {
  'https://accserv.classe.cornell.edu/svn/trunk/src/open_spacecharge' :    '50675',
  'https://accserv.classe.cornell.edu/svn/packages/PGPLOT' :               '50204',
  'https://accserv.classe.cornell.edu/svn/packages/fftw' :                 '53185',
  'https://accserv.classe.cornell.edu/svn/packages/fgsl' :                 '53176',
  'https://accserv.classe.cornell.edu/svn/packages/gnu_utilities_src' :    '51479',
  'https://accserv.classe.cornell.edu/svn/packages/gsl' :                  '53177',
  'https://accserv.classe.cornell.edu/svn/packages/hdf5' :                 '53179',
  'https://accserv.classe.cornell.edu/svn/packages/lapack' :               '53213',
  'https://accserv.classe.cornell.edu/svn/packages/lapack95' :             '53051',
  'https://accserv.classe.cornell.edu/svn/packages/openmpi' :              '51988',
  'https://accserv.classe.cornell.edu/svn/packages/plplot' :               '53180',
  'https://accserv.classe.cornell.edu/svn/packages/xraylib' :              '53181',
  'https://accserv.classe.cornell.edu/svn/trunk/util' :                    '53987',
  'https://accserv.classe.cornell.edu/svn/trunk/build_system' :            '53259',
}

all_pass = True
for repo, rev in svn_repo.items():
  p = subprocess.Popen(f"svn info {repo} | grep 'Last Changed Rev'", stdout=subprocess.PIPE, shell=True)
  (output, err) = p.communicate()
  output = output.decode("utf-8").split(':')[1].strip()
  if output == rev:
    print (f'Revision is: {output} for: {repo}')
  else:
    print (f'\033[1m\033[91mREVISION CHANGE: {output}/{rev} for: {repo}\033[0m')
    all_pass = False

if all_pass:
  print ('No SVN repo change.')
else:
  print ('\033[1m\033[91mSVN REPO CHANGE!!!!!!\033[0m')
