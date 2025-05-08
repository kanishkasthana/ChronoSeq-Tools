#!/bin/bash
#SBATCH --job-name=N701_ChronoSeq_pipeline      # Job name	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=10            # Number of CPU cores per task
#SBATCH --mem=70G                    # Job memory request
#SBATCH --time=200:00:00              # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log     # Standard output and error log
#SBATCH --exclude=compute-[4]
pwd; hostname; date
source $HOME/anaconda3/bin/activate environment_name
python /stg1/data2/kanishk/chrono-seq-tools/ChronoSeqPipelineNoCorrection.py /stg1/data2/kanishk/BaseSpace_Data/May2023QC_UsingDropSeqDevice/18156FL-30-01-01_ds.8997c1661c2c4e9aba573a1031e46a97/18156FL-30-01-01_S1_L001_R1_001.fastq.gz /stg1/data2/kanishk/BaseSpace_Data/May2023QC_UsingDropSeqDevice/18156FL-30-01-01_ds.8997c1661c2c4e9aba573a1031e46a97/18156FL-30-01-01_S1_L001_R2_001.fastq.gz /stg1/data2/kanishk/mixedHumanMouse/May2023QC_UsingDropSeqDevice/N701_DMEM_Chrono_Break N701_DMEM_Chrono_Break 

date
