#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0031   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=8.postproc        # some descriptive job name of your choice
#SBATCH --output=stdout/%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=stderr/%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=gpu         # which partition to use, default on MARS is “nodes"
#SBATCH --gres=gpu
#SBATCH --time=0-06:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=128G                # memory required per node, in the form of &#91;num]&#91;M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=32       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node


############# LOADING MODULES (optional) #############
module load apps/apptainer
module load apps/nvidia-cuda

############# MY CODE #############
apptainer run --nv helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif helixer_post_bin \
        all_wbp/h5s/heterorhabditis_bacteriophora/test_data.h5 all_wbp/heterorhabditis_bacteriophora_predictions_full.h5 \
        100 0.1 0.8 60 all_wbp/heterorhabditis_bacteriophora_helixer.gff3
