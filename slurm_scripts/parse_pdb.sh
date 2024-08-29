#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0031   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=parse_pdb        # some descriptive job name of your choice
#SBATCH --output=stdout/%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=stderr/%x-%J.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-01:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=1G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node
#SBATCH --array=1-20

############# LOADING MODULES (optional) #############
module load apps/python3

############# MY CODE #############
source .venv_wbp_scratch/bin/activate
python3 ~/wbp_scratch/parse_pdb.py folder${SLURM_ARRAY_TASK_ID}/ >> pLDDT_stats_plus.csv
deactivate
