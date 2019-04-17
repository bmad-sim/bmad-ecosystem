#!/bin/bash

# The following is for using the Cornell Compute Farm for running lux_mpi.
# For Off-Site work see the documentation for running programs on the Bmad Web site
#		https://wiki.classe.cornell.edu/ACC/ACL/RunningPrograms

# To use this script, use the command:
#   qsub run_lux_mpi.sh
# Further documentation is at:
#		https://wiki.classe.cornell.edu/Computing/GridEngine

# Set shell, but the default is /bin/bash
#$ -S /bin/bash

# To set the run time name of the job to “JOBNAME”
#$ -N JOB

# To make sure that the .e and .o file arrive in the working directory
#$ -cwd

# Set the parallel environment (pe) to our general OpenMPI parallel environment (sge_pe) and number of slots (cores) required.
# This will use multiple nodes, if the “slot” number is greater than 32.
#$ -pe sge_pe 4

# Causes Grid Engine to send a SIGUSR2/12 signal to your program - needed for proper slave job cleanup.
#$ -notify

# Enable Reservations
#$ -R y

# Enable Kerberos ticket:
#kinit -k -t /home/$USER/etc/$USER-keytab $USER

# Set OpenMPI Program
export EXEC="lux_mpi"
export BINDIR="/home/dcs16/linux_lib/production/bin"
export BINARY=${BINDIR}/${EXEC}

# These exports needed for OpenMPI
export MPIDIR=/nfs/acc/libs/Linux_x86_64_intel/devel/packages/production
export MPIBINDIR=${MPIDIR}/bin
export MPILIBDIR=${MPIDIR}/lib
export MPIRUN=${MPIBINDIR}/mpirun
export PATH=${MPIBINDIR}:${PATH}

# MPI Options for job control in SGE and using SSH as the transport.  Do Not modify.
export MPI_OPTS="-v -mca orte_forward_job_control 1 -mca orte_rsh_agent ssh"

# Set the Library Path
export LD_LIBRARY_PATH=${MPILIBDIR}:${LD_LIBRARY_PATH}

# Run the OpenMPI Executable
${MPIRUN} ${MPI_OPTS} ${BINARY}
