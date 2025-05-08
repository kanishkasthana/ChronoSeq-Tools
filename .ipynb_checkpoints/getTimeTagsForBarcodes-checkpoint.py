#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pysam
import pandas as pd
import operator
from progress.bar import Bar


# In[2]:


summaryFileName="/home/kanishk/mixedHumanMouse/star_gene_exon_tagged.dge.summary.txt"
bamFileName="/home/kanishk/mixedHumanMouse/merged_sample.sam"
dgeFileName="/home/kanishk/mixedHumanMouse/star_gene_exon_tagged.dge.txt.gz"
barcodesFileName="/home/kanishk/mixedHumanMouse/cell_Barcode_extraction_only/out_cell_readcounts.txt.gztop_barcode_list.txt"


# In[3]:


summaryFile=pd.read_table(summaryFileName,header=5)
barcodes_ordered_by_num_of_transcripts=summaryFile.CELL_BARCODE


# In[4]:


barcode_df=pd.read_table(barcodesFileName,header=None)
barcodes=barcode_df[0]


# In[5]:


bamFile=pysam.AlignmentFile(bamFileName, "r")
BamRecords=bamFile.fetch(until_eof=True)


# In[6]:


#Objects from this class will form a Key:Value pair in a dictionary with barcodes we selected
class consensus_time_tag:
    def __init__(self):
        self.firstPosition={'A':0,'T':0,'G':0,'C':0,'N':0}
        self.secondPosition={'A':0,'T':0,'G':0,'C':0,'N':0}
        self.thirdPosition={'A':0,'T':0,'G':0,'C':0,'N':0}
        
    def update(self,tag):
        #Increase count each new tag
        self.firstPosition[tag[0]]+=1
        self.secondPosition[tag[1]]+=1
        self.thirdPosition[tag[2]]+=1
    
    def get_consensus_tag(self):
        #Getting bases with most counts
        first_concensus_base=max(self.firstPosition.items(),key=operator.itemgetter(1))[0]
        second_concensus_base=max(self.secondPosition.items(),key=operator.itemgetter(1))[0]
        third_concensus_base=max(self.thirdPosition.items(),key=operator.itemgetter(1))[0]
        concensus_tag=""+first_concensus_base+second_concensus_base+third_concensus_base
        return concensus_tag


# In[7]:


#Initializing barcode dictionary. Each key value pair is Barcode: concensus_time_tag object. 
barcode_dict={}
for cell_barcode in barcodes:
    barcode_dict[cell_barcode]=consensus_time_tag()


# In[9]:


print("Counting Total Records in BAM File. Please wait.")
total_records=bamFile.count(until_eof=True)
print("There are ",total_records," records in this file.")


# In[10]:


#Query BAM file for our selected Barcodes. If the Barcodes match then we update the time-tag counts
print("Started Processing BAM file to get Time Tags for the Cell Barcodes:")
onePercentCount=int((total_records/100))
bar = Bar('Processing', max=100, suffix='%(percent)d%%')
for record in BamRecords:
    total_records+=1
    if total_records%onePercentCount==0:
        bar.next()
    print("Finished Reading ",total_records," records")
    time_tag=record.get_tag('YT')
    cell_barcode=record.get_tag('XC')
    if cell_barcode in barcode_dict:
        barcode_dict[cell_barcode].update(time_tag)             
bar.finish()


# In[11]:


#Compute the final concensus tags and make new df using summary file.
concensus_tags=[]
for cell_barcode in barcodes:
    tagBarcodePair=[cell_barcode,barcode_dict[cell_barcode].get_consensus_tag()]
    summary_vals=summaryFile.loc[summaryFile.CELL_BARCODE.str.contains(cell_barcode)].iloc[:,1:].values[0].tolist()
    concensus_tags.append(tagBarcodePair+summary_vals)


# In[12]:


concensus_tag_df=pd.DataFrame(concensus_tags,columns=["CELL_BARCODE","TIME_TAG","NUM_GENIC_READS","NUM_TRANSCRIPTS","NUM_GENES"])


# In[13]:


concensus_tag_df.to_csv(summaryFileName+".time_tags.txt",sep="\t",index=False)


# In[14]:


bamFile.close()


# In[ ]:




