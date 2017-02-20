# Script for running jobs on the grid.

#$ -S /bin/sh
#$ -l arch=linux-x64
#$ -l mem_free=6G

# Specify here where the log files go
#$ -o /home/shanksj/chess/freq_map_jobs/qlog
#$ -e /home/shanksj/chess/freq_map_jobs/qlog

#. /nfs/acc/libs/util/acc_vars.sh

. /home/shanksj/.bashrc

echo "---Job Information---"
echo "Directory:  " $workingdir
echo "Input file: " $inputfile
echo "Hostname:   " $HOSTNAME
echo "OS:         " `cat /etc/redhat-release`
echo "-------------------------------"
echo ""
echo ""

cd $workingdir

# Change this to point to your local executable.
/home/shanksj/bmad_dist/production/bin/freq_map $inputfile 
