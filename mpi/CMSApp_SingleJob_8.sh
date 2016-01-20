#!/bin/bash
## Name of the job
#PBS -N CMSApp_MPI

## Number of nodes (in this case I require 7 nodes all the slots)
## The total number of nodes passed to mpirun: -ppn 1 allow 1 instance  mpi on each nodes
#PBS -l nodes=8:ppn=2 


#! Name of output files for std output and error;
#! if non specified defaults are <job-name>.o<job number> and <job-name>.e<job-number>
#PBS -e CMSApp_MPI.$PBS_JOBID.err
#PBS -o CMSApp_MPI.$PBS_JOBID.out
#PBS -l cput=96:00:00
#PBS -l walltime=96:00:00
#PBS -q high

#PBS -W x=\"NACCESSPOLICY:SINGLEJOB\"
#PBS -l naccesspolicy=singlejob
 
## Counts the number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

 
## Create a machine file for Myrinet
#echo $NPROCS >$PBS_JOBID.nodefile
#awk '{if ($0 in vett) print $0 " " 7; else print $0 " " 6 ; vett[$0]="x"}' $PBS_NODEFILE >>$PBS_JOBID.nodefile
#awk '{print $0 ; vett[$0]="x"}' $PBS_NODEFILE >>$PBS_JOBID.nodefile
 


export OMP_NUM_THREADS=24


	echo "                                           "
	echo "=========== START CMSApp_MPI.cc   =============="
	date
	echo "==========================================="
	echo "                                           "
mpirun  -ppn 1 ./CMSApp_MPI /home/mba/pp14TeVminbias.root
   
	echo "==========================================="
	echo "                                           "

date
