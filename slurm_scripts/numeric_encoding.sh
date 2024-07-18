#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0031   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=4.geenuff2h5        # some descriptive job name of your choice
#SBATCH --output=stdout/%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=stderr/%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=gpu         # which partition to use, default on MARS is “nodes"
#SBATCH --gres=gpu
#SBATCH --time=0-06:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=64G                # memory required per node, in the form of &#91;num]&#91;M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=8       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node
#SBATCH --array=0-45

############# LOADING MODULES (optional) #############
module load apps/apptainer
module load apps/nvidia-cuda

############# MY CODE #############
for sp in `cat batches/batch${SLURM_ARRAY_TASK_ID}`
do 
        mkdir -p all_wbp/h5s/$sp
        apptainer run --nv helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif geenuff2h5.py \
        --input-db-path for_GeenuFF_new/$sp/output/$sp.sqlite3 \
        --h5-output-path all_wbp/h5s/$sp/test_data.h5
done
