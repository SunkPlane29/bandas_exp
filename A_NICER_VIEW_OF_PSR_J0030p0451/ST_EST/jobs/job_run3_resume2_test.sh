#!/bin/bash
#SBATCH -N 5
#SBATCH --tasks-per-node=24
#SBATCH -t 01:00:00
#SBATCH -p normal
#SBATCH --constraint=haswell
#SBATCH --job-name=run3_r2

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load intel/2017b
module load python/2.7.9

cp -r $HOME/NICER_analyses/J0030_ST_EST $TMPDIR

cd $TMPDIR/J0030_ST_EST

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

srun python main_run3_resume2_test.py > out_run3_resume2_test 2> err_run3_resume2_test

cp out_run3_resume2_test err_run3_resume2_test $HOME/NICER_analyses/J0030_ST_EST/.
#end of job file
