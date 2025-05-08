#Author Kanishk Asthana kasthana@eucsd.edu
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
import os
from datetime import datetime,timedelta
import argparse
import sys

#Default Values
timeTaggedDGE1File=""
timeTaggedDGE2File=""
outputFile=""

#For writing output to the stdout in realtime
def write(*something):
    print(*something)
    sys.stdout.flush()

def parse_file(input_filename):
    if not os.path.isfile(input_filename):
        raise argparse.ArgumentTypeError("File does not exist. Please use a valid file path.")
    return(input_filename)

parser=argparse.ArgumentParser(description="Script to Combine two Time-Tagged DGEs into a Single Time Series. The Maximum Value of the Time-Tag from the First DGE is summed with the Time-Tags of the Second DGE. Useful when you want to combine and multiplex the same 12 Time-Tags for longer time series with different Illumina I7 Indices. Can also be useful for Batch Correction.")
parser.add_argument("TIME_TAGGED_DGE1",help="Path for the First Time-Tagged DGE Generated using getCombinedDGEWithTimeTags.py", type=parse_file)
parser.add_argument("TIME_TAGGED_DGE2",help="Path for the Second Time-Tagged DGE Generated using getCombinedDGEWithTimeTags.py",type=parse_file)
parser.add_argument("OUTPUT_FILENAME",help="Please enter a Valid Path for the Combined Output table. Remember to end the Filename with .csv.gz",type=str)
args=parser.parse_args()
write(args)

timeTaggedDGE1File=args.TIME_TAGGED_DGE1
timeTaggedDGE2File=args.TIME_TAGGED_DGE2
outputFile=args.OUTPUT_FILENAME

timeTaggedDGE1File=os.path.abspath(os.path.expanduser(timeTaggedDGE1File))
timeTaggedDGE2File=os.path.abspath(os.path.expanduser(timeTaggedDGE2File))


write("Time Tagged DGE 1 Path: "+timeTaggedDGE1File)
write("Time Tagged DGE 2 Path: "+timeTaggedDGE2File)

script_start_time=datetime.now()

start_time=datetime.now()

write("Started Reading Both Datasets at",start_time,".")

timeTaggedDGE1=pd.read_table(timeTaggedDGE1File,sep=",")
timeTaggedDGE2=pd.read_table(timeTaggedDGE2File,sep=",")

write("Started Combining into a Single Time-Series.")

#Adding The Maximum Value of Time-Tags from First Experiment to Second
max_value=timeTaggedDGE1["FINAL TIME TAG"].max()

#Shifting Time-Tag Up by Maximum Value from First DGE File
timeTaggedDGE2["FINAL TIME TAG"]=timeTaggedDGE2["FINAL TIME TAG"]+max_value

#Getting Common Genes in Both Datasets
commonGenes=sorted(list(set(timeTaggedDGE1.columns[1:-2]) & set(timeTaggedDGE2.columns[1:-2])))

#Deleting Genes that are not Common from Dataset
timeTaggedDGE1=timeTaggedDGE1[["CELL BARCODES"]+commonGenes+["FINAL TIME TAG"]]
timeTaggedDGE1["COUNTS"]=timeTaggedDGE1.iloc[:,1:-1].sum(axis=1) #Recalculating Counts
timeTaggedDGE1=timeTaggedDGE1[["CELL BARCODES"]+commonGenes+["COUNTS","FINAL TIME TAG"]]#Reordering Columns to Match Original Configuration.

#Deleting Genes that are not Common from Dataset
timeTaggedDGE2=timeTaggedDGE2[["CELL BARCODES"]+commonGenes+["FINAL TIME TAG"]]
timeTaggedDGE2["COUNTS"]=timeTaggedDGE2.iloc[:,1:-1].sum(axis=1) #Recalculating Counts
timeTaggedDGE2=timeTaggedDGE2[["CELL BARCODES"]+commonGenes+["COUNTS","FINAL TIME TAG"]]#Reordering Columns to Match Original Configuration.

combinedTimeTaggedDGE=pd.concat([timeTaggedDGE1,timeTaggedDGE2],ignore_index=True)
combinedTimeTaggedDGE = combinedTimeTaggedDGE.drop_duplicates(subset=["CELL BARCODES"], keep=False)

write("Final Dimensions of Combined Dataset:")
write(combinedTimeTaggedDGE.shape)

#Sorting by Final Time Tag
combinedTimeTaggedDGE=combinedTimeTaggedDGE.sort_values("FINAL TIME TAG")


#Statistics for Each Time-Tag
time_tag_counts=combinedTimeTaggedDGE.groupby("FINAL TIME TAG").count().iloc[:,0:1]
time_tag_counts.to_csv(outputFile+".time_tag_counts_summary.csv",index=False)

write("Generating Bar Plot for Average number of Unique Transcripts captured per Time-Tagged Bead.")
#Plotting Average Unique Counts Detected per Cell Barcode for Each Time Tag
sns.set_theme(style="whitegrid")
ax = sns.barplot(x="FINAL TIME TAG", y="COUNTS", data=combinedTimeTaggedDGE)
ax.set(ylabel="Mean UMI/Time-Tag-Bead")
ax.set(xlabel="Time Tag Number")
fig=ax.get_figure()
fig.savefig(outputFile+".Tags_counts_bar_plot.png",dpi=600)


#Exporting Data
write("Exporting Combined Data to Compressed CSV File.")
combinedTimeTaggedDGE.to_csv(outputFile,index=False,compression="gzip")

write("Total Execution Time:",datetime.now()-script_start_time)



