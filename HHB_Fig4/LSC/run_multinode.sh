#!/bin/bash

#SBATCH -J HHB_LSC
#SBATCH -t 168:00:00
#SBATCH -n 48
#SBATCH -N 4
#SBATCH -p astra
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -x c00[55,67]
#SBACTH --exclusive
# Enter the working directory
cd ${SLURM_SUBMIT_DIR}

echo "starting $SLURM_JOB_ID at `date` on `hostname`"

echo "tmp directory is ${TMPDIR}"

echo "$USER"

MYTMP=/tmp/${USER}/${SLURM_JOB_ID}

mpiexec -pernode /usr/bin/mkdir -vp $MYTMP #|| exit $?

#echo "Copying data over... "

mpiexec -pernode scp -v ${SLURM_SUBMIT_DIR}/dyn.x  $MYTMP/. #|| exit $?
mpiexec -pernode scp -v ${SLURM_SUBMIT_DIR}/input_mD  $MYTMP/. #|| exit $?

echo "$(pwd)"

cd $MYTMP

echo "$(pwd)"

module swap gnu8 intel
#run fortran code
time mpiexec -np 48 ./dyn.x

# Copy output files back from the temp directory to working directory
mpiexec -pernode rsync -r $MYTMP/ $SLURM_SUBMIT_DIR/ #|| exit $?

rm -rf $MYTMP

exit 0 


 
