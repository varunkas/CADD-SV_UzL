#! /bin/bash

#SBATCH --partition=shortterm
#SBATCH -c 8
#SBATCH --mem=100GB
#SBATCH --job-name=CADD-SV
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

# Load your necessary modules:
conda activate run.caddsv

# Submit the Nextflow Script:
snakemake  --use-conda --configfile config.yml -j 8 --conda-prefix $WORK/cadd_sv/
