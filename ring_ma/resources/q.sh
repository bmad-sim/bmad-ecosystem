# Script for running jobs on the grid.

#$ -S /bin/csh
#$ -l arch=lx24-x86

# Specify here where the log files go
#$ -o /nfs/acc/user/helms/ctf/jobs/qlog
#$ -e /nfs/acc/user/helms/ctf/jobs/qlog

echo "Input file       : " $inputfile
cd $workingdir

# Change this to point to your local executable.
/nfs/acc/user/helms/ctf/bin/ring_ma $inputfile
