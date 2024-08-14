#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0031   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=esmfold        # some descriptive job name of your choice
#SBATCH --output=stdout/%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=stderr/%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=gpu         # which partition to use, default on MARS is “nodes"
#SBATCH --gres=gpu
#SBATCH --time=0-12:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=128G                # memory required per node, in the form of &#91;num]&#91;M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=4       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node
#SBATCH --array=1-20

############# LOADING MODULES (optional) #############
module load apps/esm

############# MY CODE #############
for x in `find ~/wbp_scratch/data/novel_orthologue_sequences/braker3/folder${SLURM_ARRAY_TASK_ID} -type f`; do
	ln -s $x fasta/
    fn="${x##*/}"
	esm-fold -i fasta/$fn -o pdb/ --chunk-size 128 > stdout/$fn.stdout
	tail -1 stdout/$fn.stdout >> esm_pLDDTs_acc.txt
done
