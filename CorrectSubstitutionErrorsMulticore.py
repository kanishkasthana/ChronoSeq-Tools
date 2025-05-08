#Author Kanishk Asthana kasthana@eng.ucsd.edu
import pysam
from datetime import datetime,timedelta
import numpy as np
import argparse
import os
import sys
from multiprocessing import Pool
import copy

#For writing output to the stdout in realtime
def write(*something):
    print(*something)
    sys.stdout.flush()

#Default Values
bamFileName=""
outBamFileName=""
MIN_UMI_COUNT=20 #Minimum UMI count for Barcode to be considered for Collapse. Not counting UMIs for simplicity.Change this later because it breaks easily. #ToDO
SUBSTITUTION_ERROR_FREQ=0.93 #Minimum frequency of majority barcode in substitution pairs needed for collapse
max_cores=1

def parse_file(input_filename):
    if not os.path.isfile(input_filename):
        raise argparse.ArgumentTypeError("File does not exist. Please use a valid file path.")
    return(input_filename)
    
def check_umi(input_umi):
    value=int(input_umi)
    assert value>0,"Please enter positive values only for UMI Count Filter!"
    return value

parser=argparse.ArgumentParser(description="Script to detect and Correct Substitution errors in the Cell Barcodes. These errors are likely sequencing errors and not Bead Synthesis Errors")
parser.add_argument("INPUT_FILENAME",help="Aligned Merged and Gene labeled BAM file for Correcting Substitution Errors.", type=parse_file)
parser.add_argument("OUTPUT_FILENAME",help="Please enter a Valid Path for the Error Corrected output BAM file.",type=str)
parser.add_argument("NUM_CORE",help="Maximum number of Cores to be used for the Computation.",type=int)
parser.add_argument("-umi","--MIN_UMI", help="Minimum UMI count per barcode to be considered for Collapse. Default Value is 20.",type=check_umi)

args=parser.parse_args()
write(args)

if args.MIN_UMI is not None:
    MIN_UMI_COUNT=args.MIN_UMI

bamFileName=args.INPUT_FILENAME
outBamFileName=args.OUTPUT_FILENAME
max_cores=args.NUM_CORE

script_start_time=datetime.now()

class CellBarcode:
    CellBarcodesWithEnoughCounts=[]
    NumberofBarcodesCollapsed=0
    def __init__(self,barcode,umi):
        self.count=1 #Number of reads with that Barcode. When initializing you count the barcode you initialize with
        self.barcode=barcode
        self.barcodeToReturn=self #This will be updated with a different reference if this barcode is merged with another
        self.umi_dict={}
        self.umi_dict[umi]=1
        
    def increase_count(self,umi):
        self.count+=1 #Increase count if you see a barcode
        if umi in self.umi_dict:
            self.umi_dict[umi]+=1
        else:
            self.umi_dict[umi]=1
            
    def computeHasEnoughCounts(self):
        if len(self.umi_dict.keys())>MIN_UMI_COUNT:
            CellBarcode.CellBarcodesWithEnoughCounts.append(self.barcode)
            
    def combineIfSubstitution(self,CellBarcode2):
        largerBarcode=self
        smallerBarcode=CellBarcode2
        if smallerBarcode.count>largerBarcode.count:
            largerBarcode=CellBarcode2
            smallerBarcode=self
        freq=largerBarcode.count/(largerBarcode.count+smallerBarcode.count)
        if freq>SUBSTITUTION_ERROR_FREQ:
            smallerBarcode.barcodeToReturn=largerBarcode
            largerBarcode.count+=smallerBarcode.count
            CellBarcode.NumberofBarcodesCollapsed+=1

barcode_dict={}

#If you get an error here your file probably not correctly formated. Make sure you have a header.
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)

start_time=datetime.now()
prevMil=start_time
write("Started Processing BAM file at",start_time,". Getting Cell Barcodes to correct Illumina Sequencing Base Substitution Errors!")

total_records=0
for record in BamRecords:
    

    #For printing progress
    total_records+=1
    if total_records%1000000==0:
        time_taken=datetime.now()-prevMil
        write("Finished processing ",total_records,"\trecords at",datetime.now(),". Previous 1000000 Records took ",time_taken.total_seconds(),"s")
        prevMil=datetime.now()
    

    #Main Logic
    cell_barcode=record.get_tag('XC')
    umi=record.get_tag('XM')
    
    if cell_barcode in barcode_dict:
        barcode_dict[cell_barcode].increase_count(umi)
    else:
        barcode_dict[cell_barcode]=CellBarcode(cell_barcode,umi)


total_time=datetime.now()-start_time
write("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

bamFile.close()

for barcode in barcode_dict.keys():
    barcode_dict[barcode].computeHasEnoughCounts()

write(len(CellBarcode.CellBarcodesWithEnoughCounts),"Cell barcodes have enough UMIs for further processing.")

#Gives total number of elements that need to be computed. 
#Should be helpful when deciding how to divide the computation.
def getTotalElementsInUpperTriangle(vector_length):
    return int((vector_length*vector_length-vector_length)/2)

#Dividing Upper Triangle into roughly equal parts for the Different Processes in Parallel
def getEqualChunksForUpperTriangle(barcodes_list,max_cores):
    vector_length=len(barcodes_list)
    cumsum=0
    chunk=getTotalElementsInUpperTriangle(vector_length)/max_cores
    starting_i=[0]
    max_val=chunk
    ending_i=[]
    for i in range(0,vector_length-1):
        cumsum+=vector_length-1-i
        if cumsum>max_val:
            max_val=cumsum+chunk
            starting_i.append(i)
            ending_i.append(i-1)
    ending_i.append(vector_length-1)
    output_pairs=[]
    for position in range(0,len(starting_i)):
        output_pairs.append((starting_i[position],ending_i[position]))
    return output_pairs

chunks=getEqualChunksForUpperTriangle(CellBarcode.CellBarcodesWithEnoughCounts,max_cores)
num_cores=len(chunks)
write("Dividing into the Following Chunks for Parallel Processing:")
write(chunks)
write("Number of Cores that will be used: ",num_cores)

pairs_list=[]

chunk_lists=[]

write("Collapsing Barcodes at a Hamming Distance of 1 and with the Major Barcode present above frequency",SUBSTITUTION_ERROR_FREQ," at",datetime.now())

#Calculating Distances for the upper triangle of the distance matrix.
def computeChunk(start,end,barcode_list):
    def hammingDistance(barcode1,barcode2): #Stop calculating as soon as distance is more than 1. Saves a lot of time and doesn't change result.
        distance=0
        i=0
        while distance<2 and i<len(barcode1):
            if barcode1[i]!=barcode2[i]:
                distance+=1
            i+=1
        return distance

    pairs_list=[]
    vector_length=len(barcode_list)
    for i in range(start,end+1):
        for j in range(i+1,vector_length):
            firstBarcode=barcode_list[i]
            secondBarcode=barcode_list[j]
            if hammingDistance(firstBarcode,secondBarcode)==1:
                pairs_list.append((i,j))
    return pairs_list

with Pool(num_cores) as pool:
    params=[]
    for start,end in chunks:
        barcodes_list_copy=copy.deepcopy(CellBarcode.CellBarcodesWithEnoughCounts)
        params.append((start,end,barcodes_list_copy))
    results=[pool.apply_async(computeChunk,p) for p in params]
    for r in results:
        chunk_lists.append(r.get())

for sub_list in chunk_lists:
    for pair in sub_list:
        pairs_list.append(pair)

for i,j in pairs_list:
    firstBarcode=CellBarcode.CellBarcodesWithEnoughCounts[i]
    secondBarcode=CellBarcode.CellBarcodesWithEnoughCounts[j]
    barcode_dict[firstBarcode].combineIfSubstitution(barcode_dict[secondBarcode])
    
write(CellBarcode.NumberofBarcodesCollapsed,"barcodes pairs were collapsed! At ",datetime.now())

#Writing to new File with Updated Barcodes
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)
outBamFile=pysam.AlignmentFile(outBamFileName,"wb",template=bamFile)

start_time=datetime.now()
prevMil=start_time
write("Started Writing Updated BAM file at",start_time,".")

total_records=0
for record in BamRecords:    

    #For printing progress
    total_records+=1
    if total_records%1000000==0:
        time_taken=datetime.now()-prevMil
        write("Finished processing ",total_records,"\trecords at",datetime.now(),". Previous 1000000 Records took ",time_taken.total_seconds(),"s")
        prevMil=datetime.now()
    
    #Main Logic
    cell_barcode=record.get_tag('XC')
    updated_barcode=barcode_dict[cell_barcode].barcodeToReturn.barcode #Getting Updated Barcode
    record.set_tag('XC',updated_barcode)
    outBamFile.write(record)

total_time=datetime.now()-start_time
write("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

bamFile.close()
outBamFile.close()

write("Total Execution Time:",datetime.now()-script_start_time)