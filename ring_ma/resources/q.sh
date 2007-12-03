#$ -S /bin/csh
#$ -l arch=lx24-x86
#$ -o /nfs/acc/user/helms/ctf/jobs/qlog
#$ -e /nfs/acc/user/helms/ctf/jobs/qlog
echo "Input file       : " $inputfile
#cd /nfs/cesr/acc/helms/ctf/jobs/2006-07-27
cd $workingdir
/nfs/acc/user/helms/ctf/bin/ring_ma $inputfile
