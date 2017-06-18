#!/usr/bin/env python

import subprocess,sys,os,re,fileinput,argparse,math
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
#repair for one sample and not two.
tempheader = open('vcf.tmp','w')
headerrewritten = open('tmp.vcf','w')
subprocess.call(['zgrep \"#\" '  + vcffi],shell=True,stdout=tempheader)
subprocess.call(['sed -e \'s/unknown\tSample1/Sample1/\' vcf.tmp'],shell=True,stdout=headerrewritten)

#--subs--#
def chooseSample(samples):

   if samples[0]['GT'] is None:
        return [samples[1]]
   else:
        return [samples[0]]


def addAnnovar(line):
    #print line.INFO
    #return NULL
    line.INFO['IDX'] = line.POS
    newsamples = chooseSample(line.samples)
    newline = vcf.model._Record(line.CHROM,line.POS,line.ID,line.REF,line.ALT,line.QUAL,line.FILTER,line.INFO,line.FORMAT,0,samples=newsamples)#,str(newsamples))
    return newline


def addFab(line):
    line.ID = line.POS
    newsamples = chooseSample(line.samples)
    newline = vcf.model._Record(line.CHROM,line.POS,line.ID,line.REF,line.ALT,line.QUAL,line.FILTER,line.INFO,line.FORMAT,0,samples=newsamples)#,str(newsamples))
    return newline


#--main--#
headertemplate = vcf.Reader(open('tmp.vcf','r'))
v = re.match('(\w+.*?)\.vcf.*',vcffi)
counter = 0
if annovartype is True:
    vcfoutfi = v.group(1) + '_2annovar.vcf'
    vcf_write = vcf.Writer(open(vcfoutfi,'w'),headertemplate)

    for record in vcf_full:
        vcf_write.write_record(addAnnovar(record))
        #(addAnnovar(record))



elif fabrictype is True:
    vcfoutfi = v.group(1) + '_2fabric.vcf'
    vcf_write = vcf.Writer(open(vcfoutfi,'w'),headertemplate)

    for record in vcf_full:
        vcf_write.write_record(addFab(record))

