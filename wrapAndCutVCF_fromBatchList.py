#!/usr/bin/env python

import subprocess,sys,os,re,fileinput,argparse,math
import vcf
import ntpath

parser = argparse.ArgumentParser(description="pulls the files and cuts them from sampleID "\t" batch ID")
parser.add_argument("--sample",help="List of the samples needed to cut.",required=True)
parser.add_argument("--bed",help="bed file of region to cut")
parser.add_argument("--wanted",help="some text snippet that would be greppable for")
parser.add_argument("--unwanted",help="some text snippet that would be greppable for")
#get through samples, get the list of hashes, then get the hashes, then grep the files
#and then send to cutter.
#

args = parser.parse_args()

samplefi = args.sample
bedfi = args.bed
regex = args.wanted
notregex = args.unwanted
#special load for vcffi? (vcffi,"r")


sample_list = open(samplefi,'r')
bedfile = open(bedfi,'r')
bedlines = []
for beds in bedfile:
    bedlines.append(beds)



#--subs--#


def grabPDH(templst):#for now just most recent, assuming that is the first returned.

    final_pdh = ""
    pdhlst = templst.split("\n")
    
    if len(pdhlst) > 2:#if valid results.
        final_pdh = pdhlst[1].split(",")[0]

    return final_pdh


def getFileWanted(pdh,rex):
    
    wantedfi = subprocess.Popen('/home/diegoisi/arvados-tools/ls_output_collections.py  '  + pdh + '  /home/diegoisi/keep/by_id/ s | egrep ' + rex,shell=True,stdout=subprocess.PIPE)
    return wantedfi.stdout.read()

def chooseSample(sid):
    retfile = None
    print sid
    filelist = subprocess.Popen('/home/diegoisi/arvados-tools/find_complete_instances.py '  + sid,shell=True,stdout=subprocess.PIPE)
    pdhmostrecent = grabPDH(filelist.stdout.read())
    full_vcf_fi = getFileWanted(pdhmostrecent,regex) 
    
    for filein in str(full_vcf_fi).split("\n"):
        searchfornotmatch = re.search(notregex,filein) 

        if searchfornotmatch is None and filein != "":
            retfile =  filein
            
    return retfile   
        
        

#--main--#

for samplename in sample_list:
    sampleid = re.search('Sample_(\w+)-EXT',samplename)

    if sampleid is not None:        
        final_file = chooseSample(sampleid.group(1))
        #print "here:" + str(final_file) + ":here"

        if final_file is not None: 
            finalvcf = vcf.Reader(open('/home/diegoisi/keep/by_id/' + final_file,'r'))
            outfilename = 'special-' + ntpath.basename(final_file)
            outfilename = re.sub(r'.gz','',outfilename)
            outvcf = vcf.Writer(open(outfilename,'w'),finalvcf)
            for bedregion in bedlines:
                reg = bedregion.split()
                cftr_region = finalvcf.fetch(reg[0],int(reg[1]),int(reg[2]))

                for variants in cftr_region:
                    outvcf.write_record(variants)
                    






