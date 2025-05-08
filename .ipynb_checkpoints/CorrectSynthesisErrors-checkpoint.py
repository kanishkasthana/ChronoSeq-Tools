#Author Kanishk Asthana kasthana@eng.ucsd.edu
import pysam
from datetime import datetime,timedelta
import argparse
import os
import pygtrie
import sys

#For writing output to the stdout in realtime
def write(*something):
    print(*something)
    sys.stdout.flush()

#Default Values
bamFileName=""
outBamFileName=""
MIN_UMI_COUNT=20 #Minimum UMI count for Barcode to be considered for Collapse. Not counting UMIs for simplicity.
num_cell_bases_missing=3 #Maximum number of Base synthesis errors from the end that will be checked for.
cell_barcode_length=12

def parse_file(input_filename):
    if not os.path.isfile(input_filename):
        raise argparse.ArgumentTypeError("File does not exist. Please use a valid file path.")
    return(input_filename)
    
def check_umi(input_umi):
    value=int(input_umi)
    assert value>0,"Please enter positive values only for UMI Count Filter!"
    return value

parser=argparse.ArgumentParser(description="Script to detect and Correct Synthesis errors in the Cell Barcodes. Please run CorrectSubstitutionErrors.py before running this script.")
parser.add_argument("INPUT_FILENAME",help="Aligned Merged, Gene labeled and Sub Corrected BAM file for Correcting Synthesis Errors.", type=parse_file)
parser.add_argument("OUTPUT_FILENAME",help="Please enter a Valid Path for the Error Corrected output BAM file.",type=str)
parser.add_argument("-umi","--MIN_UMI", help="Minimum UMI count per barcode to be considered for Collapse. Default Value is 20.",type=check_umi)
args=parser.parse_args()
write(args)

if args.MIN_UMI is not None:
    MIN_UMI_COUNT=args.MIN_UMI

bamFileName=args.INPUT_FILENAME
outBamFileName=args.OUTPUT_FILENAME

script_start_time=datetime.now()

class CellBarcode:
    CellBarcodesWithEnoughCounts=[]
    NumberofBarcodesCollapsed=0
    def __init__(self,barcode,umi):
        self.count=1 #Number of reads with that Barcode. When initializing you count the barcode you initialize with
        self.barcode=barcode
        self.wasCollapsed=False #This will be updated to True if this barcode gets collapsed.
        self.umi_dict={}
        self.umi_dict[umi]=1
        self.umi_shift=0 #This number will track number of bases to shift the UMI by when doing the correction 
        self.num_unique_umi=1
    
    #For Debugging    
    def show_info(self):
        write("Current Value:",self.barcode)
        write("Was Collapsed?",self.wasCollapsed)
        write("UMI Shift Value:",self.umi_shift)
                
    def increase_count(self,umi):
        self.count+=1 #Increase count if you see a barcode
        if umi in self.umi_dict:
            self.umi_dict[umi]+=1
        else:
            self.umi_dict[umi]=1
            
    def computeHasEnoughCounts(self):
        self.num_unique_umi=len(self.umi_dict.keys())
        if self.num_unique_umi>MIN_UMI_COUNT:
            CellBarcode.CellBarcodesWithEnoughCounts.append(self.barcode)
            
    def collapseBarcode(self,new_barcode):
        self.wasCollapsed=True
        self.umi_shift+=1
        self.barcode=new_barcode

barcode_dict={}

#If you get an error here your file probably not correctly formated. Make sure you have a header.
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)

start_time=datetime.now()
prevMil=start_time
write("Started Processing BAM file at",start_time,". Getting Cell Barcodes to correct Bead Synthesis Errors!")

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

#Building Barcode Trie:
write("Building Barcode Trie")
barcode_trie=pygtrie.CharTrie()
for barcode in CellBarcode.CellBarcodesWithEnoughCounts:
    barcode_trie[barcode]=barcode_dict[barcode]

#Recursively Update all objects in Nested Lists made in Barcode_Trie
def update(new_barcode,obj_list):
    for obj in obj_list:
        if type(obj) is CellBarcode:
            obj.collapseBarcode(new_barcode)
        else:
            update(new_barcode,obj)

write("Checking for Cell Barcode Synthesis Errors..")
for base_position in range(cell_barcode_length-1,cell_barcode_length-1-num_cell_bases_missing,-1):
    write("Checking for Synthesis Errors in Cell Barcode at Base Position ",base_position+1)
    barcodes_to_be_combined=[]
    current_barcodes=CellBarcode.CellBarcodesWithEnoughCounts
    
    #Looking for Barcodes with 4 neighbours present at base position
    for cell_barcode in current_barcodes:
        prefix=cell_barcode[0:base_position]
        try:
            #If 4 values in the list are associated with that prefix then pop the node and store it in a list
            subtrie_keys=barcode_trie.keys(prefix)
            subtrie=barcode_trie.items(prefix)
            if len(subtrie_keys)==4:
                barcodes_to_be_combined.append(subtrie)
                del barcode_trie[prefix:] #Popping Barcodes in Subtrie from Main Trie
        except KeyError: #Catching a key error if barcode already popped from the main Trie
            continue
    write(len(barcodes_to_be_combined),"subtries found at this position.")
    
    #Collapsing Barcodes
    for subtrie in barcodes_to_be_combined:
        barcode,obj_to_update=subtrie[0] #Getting First element for Generating new Barcode
        new_barcode=list(barcode)
        new_barcode[base_position]="N"#Replacing Base with N
        new_barcode="".join(new_barcode) #Converting back to String
        obj_list=[]
        for barcode,obj in subtrie:
            obj_list.append(obj)
        update(new_barcode,obj_list) #Recursively update all elements of obj_list
        barcode_trie[new_barcode]=obj_list #Adding combined Tag back to the Trie as Nested List

#Printing total number of barcodes collapsed
collapse_count=0
for barcode in barcode_dict.keys():
    if barcode_dict[barcode].wasCollapsed:
        collapse_count+=1
write(collapse_count,"barcodes were collapsed!")


def getNewUMI(original_cell_barcode,updated_barcode_object,original_umi):
    if updated_barcode_object.wasCollapsed:
        cell_barcode_tail=original_cell_barcode[-updated_barcode_object.umi_shift:]
        umi_head=original_umi[:-updated_barcode_object.umi_shift]
        newUMI=cell_barcode_tail+umi_head
        return newUMI
    else:
        return original_umi

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
    original_cell_barcode=record.get_tag('XC')
    original_umi=record.get_tag("XM")
    updated_barcode=barcode_dict[original_cell_barcode] #Getting Updated Barcode object
    new_umi=getNewUMI(original_cell_barcode,updated_barcode,original_umi)#Will return updated if changed, original if unchanged.
    record.set_tag('XC',updated_barcode.barcode)
    record.set_tag('XM',new_umi)
    outBamFile.write(record)

total_time=datetime.now()-start_time
write("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

bamFile.close()
outBamFile.close()

write("Total Execution Time:",datetime.now()-script_start_time)