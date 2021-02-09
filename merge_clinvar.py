#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
import csv
from subprocess import call
from collections import defaultdict
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
#header is:chr start   stop    rsid    homohgvs    homourl hethgvs heturl  wthgvs  wturl
parser = argparse.ArgumentParser(description="Merges the files to input for ANNOVAR. there is a 99.68% chance that this will have to be hacked.")
parser.add_argument("--varsum",help="Variant summary file from clinvar.",required=True)
parser.add_argument("--hgmd",help="HGMD file from download.",required=True)
parser.add_argument("--named",help="the date name prefix, for output files.",required=True)
#parser.add_argument("--combo",help="combination types",action='append',required=False)
#parser.add_argument("--ABthreshold",help="this is a value at which to trust allele ballance, default is 15",default=0.15)
#parser.add_argument("--lc",help="If no genotype can be found, then this is the list of novel or low coverage pages")
#parser.add_argument("--hap",help="vcf-allele-haplotypes output")
#parser.add_argument("--trans",help="the gene to transcript file for PGX",required=False)
#parser.add_argument("--rpt",help='the repeats formatted file, reuse argument for multiple searches.\n format is: \"chrom start stop base leader unit wt_version wt_length hgvs_primitive\"',action='append',required=False)
#parser.add_argument("--add",help="inject a line, usually used to add results from an external tool and add it to the combo search sections",action='append',required=False)

args = parser.parse_args()





#----------file-handling-defs----------#
#
#plan, is the main portion if file is going through, pull in other parts. fill in DOTS.
#I think OK.

hgmdfi = csv.DictReader(args.hgmd, delimiter=',')
varsum = csv.DictReader(args.varsum,delimiter=',')



#---------subroutines-----------------#





#---------------main------------------#












