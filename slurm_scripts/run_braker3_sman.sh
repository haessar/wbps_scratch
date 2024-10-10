#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0031   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=Sm_braker3        # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=3-12:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=128G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=48       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
module load apps/braker
module load apps/hisat2
module load apps/bedtools

############# MY CODE #############
braker.pl --threads 48 --gff3 --genome=$wbps_grant/fasta/schistosoma_mansoni.PRJEA36577.WBPS19.genomic_softmasked.fa --species=Schistosoma_mansoni \
     --GENEMARK_PATH=/opt/apps/alces/genemarketp/1.0.20230710/bin --PROTHINT_PATH=$wbps_grant/braker/ProtHint/bin \
     --prot_seq=$wbps_grant/fasta/WBPS_Protein_data.fa \
     --rnaseq_sets_ids=ERR11178275,ERR11178284,ERR11178301,ERR11178321,ERR11178369,ERR11178376,ERR11178378,ERR11178382,ERR11178393,ERR11178403,ERR11178413,ERR1328132,ERR1328134,ERR1328137,ERR1328141,ERR1328142,ERR1328150,ERR1328173,ERR1328181,ERR1328183,ERR1328215,ERR1328216,ERR1328220,ERR1328234,ERR1328243,ERR1328245,ERR1328258,ERR1674583,ERR1674590,ERR3489980,ERR3489989,ERR506075,ERR506105,ERR506115,ERR506116,ERR5727358,ERR5727382,ERR5727413,ERR5727432,ERR5727443,ERR5727463,ERR5727468,ERR5727473,ERR5727478,ERR5727493,ERR5727500,ERR5727503,ERR5727505,ERR5727509,ERR5727532,ERR5727537,ERR5727540,ERR5727542,ERR5727570,ERR5727571,ERR5727573,ERR5727577,ERR5727587,ERR5727591,ERR5727627,ERR5918514,ERR667063,ERR667064,ERR667098,ERR667100,ERR667106,ERR667114,ERR667123,ERR667139,ERR9854331,SRR14704933 \
     --rnaseq_sets_dirs=fastq/ 
