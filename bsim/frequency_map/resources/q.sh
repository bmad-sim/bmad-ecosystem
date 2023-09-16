# Script for running jobs on the grid.

#$ -S /bin/sh
#$ -l arch=linux-x64
#$ -l mem_free=6G

## Specify here where the log files go - MANDATORY!
## If not specified, will dump in your home dir!
#$ -o path_for_log_files
#$ -e path_for_log_files 


. ~/.bashrc

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
$ACC_EXE/frequency_map $inputfile 
