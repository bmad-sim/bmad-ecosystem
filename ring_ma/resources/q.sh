# Script for running jobs on the grid.

#$ -S /bin/sh
#$ -l arch=lx24-x86

# Specify here where the log files go
#$ -o /home/shanksj/CesrTA/ring_ma_jobs/qlog
#$ -e /home/shanksj/CesrTA/ring_ma_jobs/qlog

echo "Input file       : " $inputfile
cd $workingdir

# Change this to point to your local executable.
/home/shanksj/CesrTA/bin/ring_ma $inputfile
