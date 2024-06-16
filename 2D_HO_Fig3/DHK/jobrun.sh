#PBS -N MQCIVR
#PBS -l nodes=1:ppn=12,walltime=30:00:00
#PBS -S /bin/bash
#PBS -q default
echo $PBS_NODEFILE
echo `cat $PBS_NODEFILE`

export GMPICONF=nodeinfo/PBS_JOBID

cd $PBS_O_WORKDIR

time mpirun -np 12 ./dyn.x 

