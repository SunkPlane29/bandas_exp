#!/bin/bash
#SBATCH -J run1_NICERxXMM
#SBATCH -N 15
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=t.riley.phd@gmail.com
#SBATCH --mail-type=END

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module purge
module load python/2.7.14
module load intel/18.2.199
module load intelmpi/18.2
module load gsl/2.5-icc

dirname=$SLURM_JOBID
mkdir /tmpdir/$LOGNAME/$dirname
cp -r $HOME/J0740_STU /tmpdir/$LOGNAME/$dirname/.
cd /tmpdir/$LOGNAME/$dirname/J0740_STU/modules

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/MultiNest/MultiNest_v3.12_CMake/multinest/lib
export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo directory is $PWD
srun python $HOME/J0740_STU/modules/main_run1.py @config.ini --multinest --resume > out_run1 2> err_run1

mkdir $HOME/J0740_STU/modules/NICERxXMM/run1
mkdir $HOME/J0740_STU/modules/NICERxXMM/run1/init
mv out_run1 $HOME/J0740_STU/modules/NICERxXMM/run1/init
mv err_run1 $HOME/J0740_STU/modules/NICERxXMM/run1/init
mv samples $HOME/J0740_STU/modules/NICERxXMM/run1/samples
#end of job file
