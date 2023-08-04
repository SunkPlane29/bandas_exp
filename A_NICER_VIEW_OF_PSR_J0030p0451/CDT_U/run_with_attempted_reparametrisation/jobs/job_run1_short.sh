#!/bin/bash
#SBATCH -N 5
#SBATCH --tasks-per-node=32
#SBATCH -t 01:00:00
#SBATCH -p broadwell_short
#SBATCH --job-name=run1_short_resume2

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load intel/2017b
module load python/2.7.9

cp -r $HOME/NICER_analyses/J0030_CDTU_reparametrised $TMPDIR

cd $TMPDIR/J0030_CDTU_reparametrised

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

srun python main_run1_short.py > out_run1_short_resume2 2> err_run1_short_resume2

cp run1* out_run1_short_resume2 err_run1_short_resume2 $HOME/NICER_analyses/J0030_CDTU_reparametrised/.
#end of job file
