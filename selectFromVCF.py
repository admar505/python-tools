#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
import csv
import random
parser = argparse.ArgumentParser(description="select infor from VCF")
parser.add_argument("--vcf",help="vcf file, bgzipped and tabixed as -p vcf",required=True)
parser.add_argument("--tag",help="tag to pull",required=True,action='append')

args = parser.parse_args()
vcffi = args.vcf
taglst = args.tag

newres = open('NEW.' + vcffi + '.'  + ".txt","w")
vcf_full = vcf.Reader(open(vcffi,'r'))

#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def anyNone(rets):
    
    size = len(rets)
    not_none_ct = 0
    
    for value in rets:
        tagvals = re.split('=',value)
        if tagvals[1] != "[None]":
            not_none_ct = not_none_ct + 1

    final_return = None

    if not_none_ct <= size and not_none_ct != 0:##Whats the deal here, its flipped somehow.
        final_return = rets
        

    return final_return


def getTags(tags, varset):
    retval = []

    for tagval in tags:
        if varset[tagval] is not None:
            retval.append(tagval +  '=' + str(varset[tagval]))
    
    #make a loop or def() that checks of at least one is not none.
    return_final = anyNone(retval)
    return return_final

#####----------------MAIN--------------####      #####----------------MAIN--------------####

for line in vcf_full:
    good_tags = getTags(taglst,line.INFO)

    if good_tags is not None:#might need to manage each tag individually
        print "\t".join(good_tags)
        print_line = line.CHROM + "\t" + line.POS +  "\t" line + "\t" + "\t".join(good_tags) + "\n"
        
            





















