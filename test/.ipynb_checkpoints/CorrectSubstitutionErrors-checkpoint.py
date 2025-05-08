#Author Kanishk Asthana kasthana@eng.ucsd.edu
import pysam
from datetime import datetime,timedelta
import numpy as np
from scipy.spatial.distance import pdist
import arg
import argparse
import os

#Default Values
bamFileName=""
outBamFileName=""
MIN_READ_COUNT=100 #Minimum read count for Barcode to be considered for Collapse. Not counting UMIs for simplicity.
SUBSTITUTION_ERROR_FREQ=0.85 #Minimum frequency of majority barcode in substitution pairs needed for collapse

def parse_file(input_filename):
    if not os.path.isfile(input_filename):
        raise argparse.ArgumentTypeError("File does not exist. Please use a valid file path.")
    return(input_filename)
    
def check_rc(input_rc):
    value=int(input_rc)
    assert value>0,"Please enter positive values only for Read Count Filter!"
    return value

parser=argparse.ArgumentParser(description="Script to detect and Correct Substitution errors in the Cell Barcodes. These errors are likely sequencing errors and not Bead Synthesis Errors")
parser.add_argument("INPUT_FILENAME",help="Aligned Merged and Gene labeled BAM file for Correcting Substitution Errors.", type=parse_file)
parser.add_argument("OUTPUT_FILENAME",help="Please enter a Valid Path for the Error Corrected output BAM file.",type=str)
parser.add_argument("-rc","--MIN_RC", help="Minimum Read count per barcode to be considered for Collapse. Default Value is 100.",type=check_rc)
args=parser.parse_args()
print(args)

if args.MIN_RC is not None:
    MIN_READ_COUNT=args.MIN_RC

bamFileName=args.INPUT_FILENAME
outBamFileName=args.OUTPUT_FILENAME

script_start_time=datetime.now()

class CellBarcode:
    CellBarcodesWithEnoughCounts=[]
    NumberofBarcodesCollapsed=0
    def __init__(self,barcode):
        self.count=1 #Number of reads with that Barcode. When initializing you count the barcode you initialize with
        self.barcode=barcode
        self.barcodeToReturn=self #This will be updated with a different reference if this barcode is merged with another
    
    def increase_count(self):
        self.count+=1 #Increase count if you see a barcode
        
    def computeHasEnoughCounts(self):
        if self.count>MIN_READ_COUNT:
            CellBarcode.CellBarcodesWithEnoughCounts.append(self.barcode)
            
    def combineIfSubstitution(self,CellBarcode2):
        largerBarcode=self.barcodeToReturn
        smallerBarcode=CellBarcode2.barcodeToReturn
        if smallerBarcode.count>largerBarcode.count:
            largerBarcode=CellBarcode2.barcodeToReturn
            smallerBarcode=self.barcodeToReturn
        freq=largerBarcode.count/(largerBarcode.count+smallerBarcode.count)
        if freq>SUBSTITUTION_ERROR_FREQ:
            smallerBarcode.barcodeToReturn=largerBarcode
            largerBarcode.count+=smallerBarcode.count
            CellBarcode.NumberofBarcodesCollapsed+=1


barcode_dict={}
qname_dict={}

#If you get an error here your file probably not correctly formated. Make sure you have a header.
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)

start_time=datetime.now()
prevMil=start_time
print("Started Processing BAM file at",start_time,". Getting Cell Barcodes to correct Illumina Sequencing Base Substitution Errors!")

total_records=0
for record in BamRecords:
    

    #For printing progress
    total_records+=1
    if total_records%1000000==0:
        time_taken=datetime.now()-prevMil
        print("Finished processing ",total_records,"\trecords at",datetime.now(),". Previous 1000000 Records took ",time_taken.total_seconds(),"s")
        prevMil=datetime.now()
    

    #Main Logic
    cell_barcode=record.get_tag('XC')
    qname=record.query_name
    
    if cell_barcode in barcode_dict:
        barcode_dict[cell_barcode].increase_count()
    else:
        barcode_dict[cell_barcode]=CellBarcode(cell_barcode)
    qname_dict[qname]=barcode_dict[cell_barcode]


total_time=datetime.now()-start_time
print("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

bamFile.close()

for barcode in barcode_dict.keys():
    barcode_dict[barcode].computeHasEnoughCounts()

print(len(CellBarcode.CellBarcodesWithEnoughCounts),"Cell barcodes have enough read counts for further processing.")

def hammingDistance(barcode1,barcode2):
    distance=0
    barcode1=barcode1[0]
    barcode2=barcode2[0]
    for i in range(0,len(barcode1)):
        if barcode1[i]!=barcode2[i]:
            distance+=1
    return distance

print("Calculating Pairwise Hamming Distances between Barcodes! This may take some time. O(n^2) operation")
barcodes_array=np.array(CellBarcode.CellBarcodesWithEnoughCounts)
barcodes_array=barcodes_array.reshape(-1,1)
compressed_distances=pdist(barcodes_array,hammingDistance)

print("Collapsing Barcodes at a Hamming Distance of 1 and with the Major Barcode present above frequency",SUBSTITUTION_ERROR_FREQ)
#Getting Truth vector for all Barcode pairs with Hamming Distance of 1
truth_vector=compressed_distances==1
#Getting Indices for Barcodes pairs in the Truth Vector
indices=np.triu_indices(len(CellBarcode.CellBarcodesWithEnoughCounts),1)#Use diagnal offset=1. We don't want distance of barcodes with themself
for i in range(0,len(truth_vector)):
    if truth_vector[i]:
        firstBarcode=CellBarcode.CellBarcodesWithEnoughCounts[indices[0][i]]
        secondBarcode=CellBarcode.CellBarcodesWithEnoughCounts[indices[1][i]]
        barcode_dict[firstBarcode].combineIfSubstitution(barcode_dict[secondBarcode])

print(CellBarcode.NumberofBarcodesCollapsed,"barcodes pairs were collapsed!")

#Writing to new File with Updated Barcodes
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)
outBamFile=pysam.AlignmentFile(outBamFileName,"wb",template=bamFile)

start_time=datetime.now()
prevMil=start_time
print("Started Writing Updated BAM file at",start_time,".")

total_records=0
for record in BamRecords:    

    #For printing progress
    total_records+=1
    if total_records%1000000==0:
        time_taken=datetime.now()-prevMil
        print("Finished processing ",total_records,"\trecords at",datetime.now(),". Previous 1000000 Records took ",time_taken.total_seconds(),"s")
        prevMil=datetime.now()
    
    #Main Logic
    qname=record.query_name
    updated_barcode=qname_dict[qname].barcodeToReturn.barcode #Getting Updated Barcode for Query Name
    record.set_tag('XC',updated_barcode)
    outBamFile.write(record)

total_time=datetime.now()-start_time
print("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

bamFile.close()
outBamFile.close()

print("Total Execution Time:",datetime.now()-script_start_time)