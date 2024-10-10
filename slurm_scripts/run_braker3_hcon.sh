#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0031   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=Hc_braker3        # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=6-12:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=256G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=48       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
module load apps/braker
module load apps/hisat2
module load apps/bedtools

############# MY CODE #############
braker.pl --threads 48 --gff3 --genome=$wbps_grant/fasta/haemonchus_contortus.PRJEB506.WBPS19.genomic_softmasked.fa --species=Haemonchus_contortus \
     --GENEMARK_PATH=/opt/apps/alces/genemarketp/1.0.20230710/bin --PROTHINT_PATH=$wbps_grant/braker/ProtHint/bin \
     --prot_seq=$wbps_grant/fasta/WBPS_Protein_data.fa \
     --rnaseq_sets_ids=ERR411602,ERR411603,ERR411604,ERR411605,ERR411606,ERR411607,ERR411608,ERR411609,ERR411610,ERR411611,ERR411612,ERR411613,ERR411614,ERR411615,ERR411616,ERR411617,ERR411618,ERR411619,SRR12638992,SRR12638993,SRR12638994,SRR12638995,SRR12638996,SRR12638997,SRR12638998,SRR12638999,SRR12639000,SRR12639001,SRR12639002,SRR12639003,SRR14054472,SRR14054473,SRR14054474,SRR14054475,SRR14054476,SRR14054477,SRR14054478,SRR14054479,SRR14054480,SRR14054481,SRR14054483,SRR14054484,SRR14054485,SRR14054486,SRR14054487,SRR14054488,SRR14054489,SRR14054490,SRR14054491,SRR14262044,SRR14262045,SRR14262046,SRR14262047,SRR14262048,SRR14262049,SRR14262050,SRR14262051,SRR14262052,SRR14262053,SRR14262054,SRR14262055,SRR1616954,SRR1616956,SRR16757895,SRR16757896,SRR16757897,SRR16757898,SRR16757899,SRR16757900,SRR16757901,SRR16757902,SRR16757903,SRR21617129,SRR21617130,SRR21617131,SRR23974015,SRR23974016,SRR24491039,SRR24491040,SRR24491041,SRR24491042,SRR24491043,SRR24491044,SRR27211020,SRR27211021,SRR6325479,SRR6865584,SRR6865585,SRR6865586,SRR6865587,SRR6865591,SRR6865602,SRR6865603,SRR6865604,SRR6865605,SRR6865606,SRR6865607,SRR6865608,SRR6865609,SRR6865610,SRR6865611,SRR928055,SRR928056,SRR928057,SRR928058,SRR928059,SRR928060,SRR928061,SRR928062,SRR928063 \
     --rnaseq_sets_dirs=fastq/ 
