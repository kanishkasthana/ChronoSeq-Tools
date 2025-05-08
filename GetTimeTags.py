#Author Kanishk Asthana kasthana@eng.ucsd.edu
import pysam
from datetime import datetime,timedelta
import argparse
import os
import operator
import pandas as pd
import sys
import re

#For writing output to the stdout in realtime
def write(*something):
    print(*something)
    sys.stdout.flush()

#Default Values
bamFileName=""
outFileName=""
MIN_TT_COUNT=20 #Minimum number of Times you see a Time-Tag for a particular Cell Barcode.

def parse_file(input_filename):
    if not os.path.isfile(input_filename):
        raise argparse.ArgumentTypeError("File does not exist. Please use a valid file path.")
    return(input_filename)
    
def check_ttc(input_tt):
    value=int(input_tt)
    assert value>0,"Please enter positive values only for Time Tag Count Filter!"
    return value

parser=argparse.ArgumentParser(description="Script to get the Corrected Time-Tags for Cell Barcodes. Please run CorrectSubstitutionErrors.py and CorrectSynthesisErrors.py before running this script.")
parser.add_argument("INPUT_FILENAME",help="Aligned Merged, Gene labeled ,Sub and Syn Error Corrected BAM file", type=parse_file)
parser.add_argument("OUTPUT_FILENAME",help="Please enter a Valid Path for the output table with the Barcodes and Time-Tags.",type=str)
parser.add_argument("-ttc","--MIN_TTC", help="Minimum Time Tag count detected per barcode to be considered for Time Tag calculation. Default Value is 20.",type=check_ttc)
args=parser.parse_args()
write(args)

if args.MIN_TTC is not None:
    MIN_TT_COUNT=args.MIN_TTC

bamFileName=args.INPUT_FILENAME
outFileName=args.OUTPUT_FILENAME

script_start_time=datetime.now()

#source https://www.idtdna.com/pages/products/custom-dna-rna/mixed-bases
letter_dict={'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W', 'H': 'D', 'B': 'V', 'V': 'B', 'D': 'H', 'N': 'N', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
regex_dict={'R':"[AG]",'Y':"[CT]",'M':"[AC]",'K':"[GT]",'S':"[GC]",'W':"[AT]",'H':"[ACT]",'B':"[GCT]",'V':"[ACG]",'D':"[AGT]",'N':"[ACTG]",'A':'A','G':'G','C':'C','T':'T'}
#Time-Tag Sequences Synthesized by IDT and the Name you gave them.

time_tag_oligos={
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATNTHGHGBAAAAAAAAAAAAAAAAAAA": "SEQ1_TTGG",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACDCDTNTBAAAAAAAAAAAAAAAAAAA": "SEQ2_CCTT",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGHGHANABAAAAAAAAAAAAAAAAAAA": "SEQ3_GGAA",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATNTDCDCBAAAAAAAAAAAAAAAAAAA": "SEQ4_TTCC",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATNTSANABAAAAAAAAAAAAAAAAAAA": "SEQ5_TTAA",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATVTVTVTBAAAAAAAAAAAAAAAAAAA": "SEQ6_TTTT",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACDCKANABAAAAAAAAAAAAAAAAAAA": "SEQ7_CCAA",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACDCWGHGBAAAAAAAAAAAAAAAAAAA": "SEQ8_CCGG",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACDCDCDCBAAAAAAAAAAAAAAAAAAA": "SEQ9_CCCC",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGHGMTNTBAAAAAAAAAAAAAAAAAAA": "SEQ10_GGTT",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGHGWCDCBAAAAAAAAAAAAAAAAAAA": "SEQ11_GGCC",
 "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGHGHGHGBAAAAAAAAAAAAAAAAAAA": "SEQ12_GGGG"
}

write("Time-Tag Oligos Synthesized by IDT:")
write(time_tag_oligos)
def getReverseCompliment(inputString):
    outputString=""
    for letter in inputString[::-1]:
        outputString+=letter_dict[letter]
    return(outputString)
    
def getRegexSearchString(inputString):
    outputString=""
    for letter in inputString:
        outputString+=regex_dict[letter]
    return(outputString)

def getRegexSearchDictForOligos(oligo_dict):
    regexOligoSearchDict={}
    for key in oligo_dict:
        condensed_key="AAAAAA"+key.strip('A')+"AAAAAA"
        regexOligoSearchDict[getRegexSearchString(getReverseCompliment(condensed_key))]=oligo_dict[key]
    return(regexOligoSearchDict)

regexOligoSearchDict=getRegexSearchDictForOligos(time_tag_oligos)
write("We Will be Searching for the Following Strings for Time-Tags:")
write(regexOligoSearchDict)

def convertRegexDictToTimeTagDict(regex_dict):
    new_dict={}
    for regex in regex_dict:
        new_dict[regexOligoSearchDict[regex]]=regex_dict[regex]
    return(new_dict)

class CellBarcode:
    CellBarcodesWithEnoughCounts=[]
    def __init__(self,barcode,time_tag):
        self.barcode=barcode
        self.total_time_tag_count=0
        self.time_tag_counts_dict={}
        #Initializing Counts for Each Time Tag for a Specific Cell Barcode.
        for key in regexOligoSearchDict:
            self.time_tag_counts_dict[key]=0
        self.update(time_tag)
         
    def update(self,time_tag):
        result=None
        for key in regexOligoSearchDict:
            result=re.search(key,time_tag)
            if result is not None:
                self.total_time_tag_count+=1
                self.time_tag_counts_dict[key]+=1
                
    def computeHasEnoughCounts(self):
        if self.total_time_tag_count>MIN_TT_COUNT:
            CellBarcode.CellBarcodesWithEnoughCounts.append(self.barcode)

    def getFinalTimeTag(self):
        tag_with_max_counts=max(self.time_tag_counts_dict,key=self.time_tag_counts_dict.get)
        percentage_detected=self.time_tag_counts_dict[tag_with_max_counts]/self.total_time_tag_count
        if percentage_detected>=0.70:
            self.final_time_tag=regexOligoSearchDict[tag_with_max_counts]
        else:
            self.final_time_tag=regexOligoSearchDict[tag_with_max_counts]+" . Time Tag Collision Detected. Primary Tag is present less than 70% of total Time Tags detected."

barcode_dict={}

#If you get an error here your file probably not correctly formated. Make sure you have a header.
bamFile=pysam.AlignmentFile(bamFileName,"rb")
BamRecords=bamFile.fetch(until_eof=True)

start_time=datetime.now()
prevMil=start_time
write("Started Processing BAM file at",start_time,".")

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
    time_tag=record.get_tag('YT')
    
    if cell_barcode in barcode_dict:
        barcode_dict[cell_barcode].update(time_tag)
    else:
        barcode_dict[cell_barcode]=CellBarcode(cell_barcode,time_tag)

total_time=datetime.now()-start_time
write("Finished processing BAM file at ",datetime.now(),". Total time taken ",total_time)

bamFile.close()

for barcode in barcode_dict:
    barcode_dict[barcode].computeHasEnoughCounts()

write(len(CellBarcode.CellBarcodesWithEnoughCounts),"Cell barcodes have enough Time Tags Detected for further processing.")

barcodes_dict_df={
 "CELL BARCODES":[],
 "TIME TAG COUNTS":[],
 "FINAL TIME TAG":[]
}

for oligo in time_tag_oligos:
    barcodes_dict_df[time_tag_oligos[oligo]]=[]
write(barcodes_dict_df)

for barcode in CellBarcode.CellBarcodesWithEnoughCounts:
    cell_barcode=barcode_dict[barcode]
    cell_barcode.getFinalTimeTag()
    barcodes_dict_df["CELL BARCODES"].append(barcode)
    barcodes_dict_df["TIME TAG COUNTS"].append(cell_barcode.total_time_tag_count)
    barcodes_dict_df["FINAL TIME TAG"].append(cell_barcode.final_time_tag)
    time_tag_dict=convertRegexDictToTimeTagDict(cell_barcode.time_tag_counts_dict)
    for time_tag in time_tag_dict:
        barcodes_dict_df[time_tag].append(time_tag_dict[time_tag])

new_df=pd.DataFrame.from_dict(barcodes_dict_df)

new_df.to_csv(outFileName,index=False,sep="\t")

write("Total Execution Time:",datetime.now()-script_start_time)