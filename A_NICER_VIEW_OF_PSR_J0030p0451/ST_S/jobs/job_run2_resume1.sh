#!/bin/bash
#SBATCH -N 5
#SBATCH --tasks-per-node=32
#SBATCH -t 05:00:00
#SBATCH -p broadwell
#SBATCH --job-name=run2_r1

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load intel/2017b
module load python/2.7.9

cp -r $HOME/NICER_analyses/J0030_STS $TMPDIR

cd $TMPDIR/J0030_STS

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

srun python main_run2_resume1.py > out_run2_resume1 2> err_run2_resume1

cp run2* out_run2_resume1 err_run2_resume1 $HOME/NICER_analyses/J0030_STS/.
#end of job file