#!/bin/bash
#$ -S /bin/bash         #-- the shell for the job
#$ -q gpu.q             #-- use the gpu queue
#$ -j y                 #-- tell the system that the STDERR and STDOUT should be joined
#$ -cwd                 #-- tell the job that it should start in your working directory
#$ -l mem_free=2G       #-- submits on nodes with enough free memory
#$ -l h_rt=3:00:00      #-- runtime limit - max 2 weeks == 336 hours
##$ -l hostname=!('qb3-atgpu*'|'qb3-atgpu**'|'qb3-iogpu4')
#$ -N eq_cftr
##$ -pe mpi_onehost 4     #we want 4 mpi nodes on 1 gpu
#$ -t 1-20              # number to run at a time
#$ -tc 1              # how many to run at a time

date
hostname

## GPU/CPU stuff
GMX_CPUS=3
export CUDA_VISIBLE_DEVICES=$SGE_GPU
export OMP_NUM_THREADS=8
echo $CUDA_VISIBLE_DEVICES
echo $OMP_NUM_THREADS

# ## GPU/CPU stuff
# GMX_CPUS=8
# export CUDA_VISIBLE_DEVICES=$SGE_GPU
# export OMP_NUM_THREADS=$GMX_CPUS
# export GMX_GPU_DD_COMMS=true
# export GMX_GPU_PME_PP_COMMS=true
# export GMX_FORCE_UPDATE_DEFAULT_GPU=true
# echo $CUDA_VISIBLE_DEVICES
# echo $OMP_NUM_THREADS


# ## Get software
# module load mpi
# module load Sali
# module load cuda/10.0.130

GMX_VER=2020.6
GMX_BIN=/wynton/home/grabe/shared/gromacs/gromacs-${GMX_VER}_CUDA10_SSE4/bin/gmx


#for full gpu offloading
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

upperdir="/wynton/home/grabe/csheen/cftr-project/cftr-10c-glpg-unbinding/run01"
 
readarray -t run_array < ${upperdir}/run_scripts/run_commands.txt
#run_array=( $(ls -t run_scripts/) )
echo $run_array
num_reps=${#run_array[@]}
runscript=${upperdir}/run_scripts/${run_array[$((SGE_TASK_ID-1))]}

#check we have the right number of jobs
#add later lol


#runscript=$1
echo $runscript

source $runscript

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
