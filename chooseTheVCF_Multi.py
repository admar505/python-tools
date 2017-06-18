#!/usr/bin/env python

import sys,os,re,fileinput,argparse,math
import vcf

parser = argparse.ArgumentParser(description="adds the position info into either the ID col or in INFO")
parser.add_argument("--vcf",help="the vcf file, has to be bgzipped and indexed",required=True)
parser.add_argument("--annovar",help="prep the genome for use in annovar workflow",action='store_true',default=False)
parser.add_argument("--fabric",help="prep the vcf for use in Fabric Workflow",action='store_true',default=False)

args = parser.parse_args()

vcffi = args.vcf
annovartype = args.annovar
fabrictype = args.fabric
#special load for vcffi? (vcffi,"r")

vcf_full = vcf.Reader(open(vcffi,'r'))


#hom_total = 0;

#--subs--#
def chooseSample(samples):

    if samples[0]['GT'] == ".":

        return (samples[1]['GT'],samples[1]['DP'],samples[1]['GT'])
    else:
        
        return (samples[0]['GT'],samples[0]['DP'],samples[0]['GT'])
    

def addAnnovar(line):
    #print line.INFO
    #return NULL
    line.INFO['IDX'] = line.POS
    
    if line.samples[0]['GT'] == ".":

        #newsamples = vcf.model._Call(line.samples[1]['GT'],line.samples[1]['DP'],line.samples[1]['GT'])  
        #newline = vcf.model._Record(line.CHROM,line.POS,line.ID,line.REF,line.ALT,line.QUAL,line.FILTER,line.INFO,line.FORMAT,0,samples=newsamples)#,str(newsamples))
        #return newline
        return "NULL"
    else:
        #newsamples = vcf.model._Call(line.samples[0]['GT'],line.samples[0]['DP'],line.samples[0]['GT'])
        #newline = vcf.model._Record(line.CHROM,line.POS,line.ID,line.REF,line.ALT,line.QUAL,line.FILTER,line.INFO,line.FORMAT,0,samples=newsamples)#,str(newsamples))
        #return newline
        return "NULL"


def addFab(line):
    line.ID = line.POS
    line.samples[0] = (line.samples)
    return line



#--main--#
v = re.match('(\w+.*?)\.vcf.*',vcffi)
counter = 0
print vcf_full.metadata
if annovartype is True:
    vcfoutfi = v.group(1) + '_2annovar.vcf'
    vcf_write = vcf.Writer(open(vcfoutfi,'w'),vcf_full)
    
   # for record in vcf_full:
    #    vcf_write.write_record(addAnnovar(record))
        
    
    
elif fabrictype is True:
    vcfoutfi = v.group(1) + '_2fabric.vcf'
    vcf_write = vcf.Writer(open(vcfoutfi,'w'),vcf_full)
   # for record in vcf_full:
   #     vcf_write.write_record(addFab(record))





