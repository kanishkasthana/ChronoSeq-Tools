#!/bin/bash
#SBATCH --job-name=STAR      # Job name	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=128G                    # Job memory request
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
pwd; hostname; date
source $HOME/anaconda3/bin/activate environment_name
STAR   --runMode genomeGenerate   --runThreadN 12   --genomeDir /stg1/data2/kanishk/mm10/starGenomeIndex   --genomeFastaFiles /stg1/data2/kanishk/mm10/mm10.fasta --sjdbGTFfile /stg1/data2/kanishk/mm10/mm10.gtf   --sjdbOverhang 149
###
date
