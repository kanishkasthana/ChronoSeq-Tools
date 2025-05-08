import argparse
import os

MIN_READ_COUNT=20

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

print(bamFileName,outBamFileName,MIN_READ_COUNT)