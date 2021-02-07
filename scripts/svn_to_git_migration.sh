# Below is the script that was used to migrate tao/python in the SVN repository to GitHub.
# Get users file:
svn log --xml --quiet | grep author | sort -u | perl -pe 's/.*>(.*?)<.*/$1 = /' > users.txt
# Manually edit the users.txt to add name and email.
vim users.txt
# Invoke the migration command
git svn clone https://accserv.lepp.cornell.edu/svn/trunk/src/tao/python --username=hslepicka --authors-file=users.txt --no-metadata --prefix "" -s pytao-git
cd pytao-git
git remote add upstream https://github.com/bmad-sim/pytao.git
git push origin master