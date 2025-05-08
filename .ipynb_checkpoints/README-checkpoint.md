# ChronoSeq Tools

These scripts are used in conjuction with Dropseq tools and picard to get the Time tags and  Digital Gene Expression Matrices for ChronoSeq data. Please follow the instructions below to get started:

+ Please download:
    + And unzip [Drop-Seq Tools](https://github.com/broadinstitute/Drop-seq/releases/download/v2.4.0/Drop-seq_tools-2.4.0.zip)  
    + The [Picard](https://github.com/broadinstitute/picard/releases/download/2.23.9/picard.jar) jar file.

+ Download and install [Anaconda](https://www.anaconda.com/download) for your computing environment.
+ The environment yml files can be found in the [config_files](config_files) folder. Please create an environment identical to [chronoseqTools.env.yml](config_files/chronoseqTools.env.yml) [using the file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).
+ [Activate the environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment) you just created. Now you can run the scripts in this repo.
> Optional Note: [This page](https://docs.icer.msu.edu/Using_conda/) or [this page](https://bioinformatics.uef.fi/guides/slurm/anaconda) might be helpful if you are running conda on [SLURM](https://slurm.schedmd.com/documentation.html).
+ Every script has a help prompt. You can get the help prompt by running:``` python <script_name> -h ``` 
> Note: Please change the settings and paths at the top of each script before you run the script. You can find these settings above a long line of ###### symbols before the import statements.
+ You will need to download the prepared alignment data for Drop-seq before you can run the pipeline from https://mccarrolllab.org/dropseq/ under "Data Resources".
> Read the [Drop-seq alignment cookbook](https://github.com/broadinstitute/Drop-seq/blob/master/doc/Drop-seq_Alignment_Cookbook.pdf) for more information.
+ You will also need to generate your Star Genome Index before you get started. The [useful_commands.txt](useful_commands.txt) or [this example SBATCH script](ExampleScriptSBATCHStarGenerateGenomeIndex.sh) can help you get started with this.
> OPTIONAL:  You can [download and Setup Basespace](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview) for your Computing Environment. This can make it easier to get Sequencing Data if it is directly uploaded to Basespace. See the 
[Useful Commands File](useful_commands.txt) for more details.

+ The compute cluster we developed these scripts for did not have SSDs for the home directories. We only had SSDs for the scratch space on each compute node. As a result the main scripts have a python [Remote Locks](https://docs.python.org/3/library/threading.html#lock-objects) enabled by default. This has several advantages:
    + This prevents too many I/O operations on the home directory storage drives which can slow down the cluster.
    + You can copy and run the analysis quicker on the SSDs.
    + You prevent multiple scripts and batch jobs from copying the same or a lot of files at the same time. This can also prevent the cluster from slowing down. 
    + If you don't want this behavior enabled then just set the ```copy_to_scratch``` variable to ```False```.
    + You will need to run the [StartLock.py](StartLock.py) script on one of the Compute Nodes to enable this feature. This way all the scripts on the batch job know they have to talk to this Compute Node to release or capture the Lock. You also need to specify this Compute node in your script by changing the value in the ```ip_or_hostname_for_remote_lock``` variable.
+ Once everything above has been sorted, you can generate the Digital Gene Expression Matrices for your sequencing data by running the following Scripts in order:
    1. [ChronoSeqPipelineNoCorrectionTimeTag.py](ChronoSeqPipelineNoCorrection.py)
    2. [getTopBarcodes.py](getTopBarcodes.py)
    3. [GetDGE.py](GetDGE.py) or [MakeMixedSpeciesPlot.py](MakeMixedSpeciesPlot.py)
> Run these Scripts inside a SBATCH script if using SLURM. Also remember to activate your environment if its different from your base environment. This example [SBATCH Script](ExampleScriptChronoSeqPipeline_SBATCH_SLURM.sh) should help.
+ [SpeciesMixingAnalysisExample.ipynb](SpeciesMixingAnalysisExample.ipynb) can help generate your Mixed Species Baryard Plots.
