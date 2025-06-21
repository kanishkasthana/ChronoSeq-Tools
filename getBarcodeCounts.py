#Author Kanishk Asthana kasthana@eng.ucsd.edu
import pysam
from datetime import datetime,timedelta
import argparse
import os
import pandas as pd
import sys

#For writing output to the stdout in realtime
def write(*something):
    print(*something)
    sys.stdout.flush()

#Default Values
bamFileName=""
outFileName=""

def parse_file(input_filename):
    if not os.path.isfile(input_filename):
        raise argparse.ArgumentTypeError("File does not exist. Please use a valid file path.")
    return(input_filename)

parser=argparse.ArgumentParser(description="Script to get a list of all Cell Barcodes and their counts in the BAM file. Barcodes are listed in descending order.")
parser.add_argument("INPUT_FILENAME",help="BAM file with merged tags after alignment", type=parse_file)
parser.add_argument("OUTPUT_FILENAME",help="Filename for Compressed output. Please suffix filename with .gz to avoid confusion.",type=str)
args=parser.parse_args()
write(args)

bamFileName=args.INPUT_FILENAME
outFileName=args.OUTPUT_FILENAME

script_start_time=datetime.now()

barcode_dict={}

#If you get an error here your file probably not correctly formated. Make sure you have a header.
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)

start_time=datetime.now()
prevMil=start_time
write("Started Processing BAM file at",start_time,". Getting Cell Barcode counts.")

total_records=0
for record in BamRecords:
    

    #For printing progress
    total_records+=1
    if total_records%1000000==0:
        time_taken=datetime.now()-prevMil
        write("Finished processing ",total_records,"\trecords at",datetime.now(),". Previous 1000000 Records took",time_taken.total_seconds(),"s")
        prevMil=datetime.now()
    

    #Main Logic
    cell_barcode=record.get_tag('XC')
    
    if cell_barcode in barcode_dict:
        barcode_dict[cell_barcode]+=1
    else:
        barcode_dict[cell_barcode]=1


total_time=datetime.now()-start_time
write("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

barcode_df=pd.DataFrame(list(barcode_dict.items()),columns=["CELL_BARCODE","COUNTS"])
barcode_df=barcode_df.sort_values("COUNTS",ascending=False).reset_index(drop=True)
barcode_df.to_csv(outFileName,index=False,header=False,sep="\t",compression='gzip')

write("Total Execution Time:",datetime.now()-script_start_time)
