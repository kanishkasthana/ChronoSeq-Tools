#!/bin/bash
#SBATCH --job-name=star_alignment      # Job name	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --mem=128G                    # Job memory request
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
pwd; hostname; date

STAR --runThreadN 12 --genomeDir /home/kanishk/Temp/starGenomeIndex --readFilesIn /home/kanishk/mixedHumanMouse/Media_Optimization_Library2/Media_Optimization_Library2_unaligned.tagged.Cell.Molecular.Time.Filtered.Adapter_Trimmed.PolyA_Trimmed.fastq --outFileNamePrefix star


date
