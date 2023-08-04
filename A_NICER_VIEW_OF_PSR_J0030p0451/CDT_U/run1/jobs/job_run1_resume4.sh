#!/bin/bash
#SBATCH -N 5
#SBATCH --tasks-per-node=32
#SBATCH -t 5-00:00:00
#SBATCH -p broadwell
#SBATCH --job-name=run1_resume4

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load intel/2017b
module load python/2.7.9

cp -r $HOME/NICER_analyses/J0030_CDTU $TMPDIR

cd $TMPDIR/J0030_CDTU

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

srun python main_run1_resume4.py > out_run1_resume4 2> err_run1_resume4

cp run1* out_run1_resume4 err_run1_resume4 $HOME/NICER_analyses/J0030_CDTU/.
#end of job file