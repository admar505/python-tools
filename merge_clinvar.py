#!/usr/bin/env python3
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
csv.field_size_limit(sys.maxsize)



#all keys will be for the 
#KEY for anything is chr-start-stop-ref-alt


#----------file-handling-defs----------#
#
#plan, is the main portion if file is going through, pull in other parts. fill in DOTS.
#I think OK.



hgmdfi = csv.DictReader(open(args.hgmd,'r'), delimiter=',')#keeps file handle permanently in memory, I dont think this is a big deal.
hgmd = {} #store HGMD lines here.

varsumfi = csv.DictReader(open(args.varsum,'r'),delimiter='\t')

varsum_out_name = "hg19_clinvar_" + args.named + ".txt"
varsum_final = open(varsum_out_name,"w")



#---------subroutines-----------------#

def getCLN(varline):#there is a double key, one column has to be split.
                        #this is almost the sole reason there is a key.
    
    varlst = []
    varlst = varline['PhenotypeIDS'].split(",")
    CLNDSDB = {}
    CLNDSDBID = {}

    for clnln in varlst:
        
        array = clnln.split(":")#because clinvar SUCKS, the ID is only the last entry in this record.
        arrct  = 0

        while arrct < (len(array) - 1):
        
            CLNDSDB[array[arrct]] = arrct
            arrct = arrct + 1


        CLNDSDBID[array[-1]] = array[-1]
    
    returnstring = ",".join(CLNDSDB.keys()) + "\t" + ",".join(CLNDSDBID.keys())
    
    return(returnstring)

def checkGuidelines(glines):#simply checks if guidelines are dashes, it needs to be dots.
    
    if "-" in str(glines):
        glines = "."

    return(glines)

def fixStrand(hgmdln):

    def __revit__(base,pos):
        baseret =  base.translate(str.maketrans({'A':'T','G':'C','T':'A','C':'G'}))
        return(baseret)
        

    if "-" in str(hgmdln['gene_strand']):
        keyh = hgmdln['chromosome'] + "-" + hgmdln['hg19_start'] + "-" + hgmdln['hg19_end'] + "-" + __revit__(hgmdln['hgmd_ref'],hgmdln['hg19_start']) + "-" + __revit__(hgmdln['hgmd_alt'],hgmdln['hg19_start'])

    else:
        keyh = hgmdln['chromosome'] + "-" + hgmdln['hg19_start'] + "-" + hgmdln['hg19_end'] + "-" + hgmdln['hgmd_ref'] + "-" + hgmdln['hgmd_alt'] 

    return(keyh) 
    



#---------------main------------------#

for hln in hgmdfi:#HGMD is going in.

    key = fixStrand(hln)
    hgmd[key] = hln['hgmd_accession'] + "\t"  + hln['variant_class'] + "\t"  + hln['primary_pubmed']

varsum_final.write("#Chr\tStart\tEnd\tRef\tAlt\tCLNSIG\tCLNDBN\tCLNACC\tCLNDSDB\tCLNDSDBID\tReview\tGuidelines\tVariationID\thgmdID\thgmdVC\thgmdPUBMED\n")


"""
for testkey in hgmd:

    print( testkey )


"""
removehgmd = {}#this is to remove the data after the sew on to the hgmd. it has to be done in phases, so that 
               #the keys can be used multiple times.


for var in varsumfi:#more about the VARSUM file.

    retrievalkey = var['Chromosome'] +"-"+ var['Start'] + "-" + var['Stop'] + "-" + var['ReferenceAlleleVCF'] + "-" + var['AlternateAlleleVCF']
    
    clnDBS = getCLN(var)

    try:
        hgmdinfo = hgmd[retrievalkey]
        removehgmd[retrievalkey] = retrievalkey # delete these records prior to print.

    except KeyError as e:
        hgmdinfo = ".\t.\t."


    changeguidelinetodot = checkGuidelines(var['Guidelines'])

    varsum_final.write(var['Chromosome'] +"\t"+ var['Start'] + "\t" + var['Stop'] + "\t" + var['ReferenceAlleleVCF'] + "\t" + var['AlternateAlleleVCF']\
             + "\t" + var['ClinicalSignificance'] + "\t" + var['PhenotypeList'] + "\t" + var['RCVaccession'] + "\t" + clnDBS + "\t" + var['ReviewStatus']\
             + "\t" + changeguidelinetodot + "\t" + var['#AlleleID'] + "\t" + hgmdinfo  +  "\n")



##Now, fill in the items which did not have a match in the clinvar files.

for cleanup in hgmd:

    if cleanup not in removehgmd:
        


        chrinf = cleanup.split("-")
        ind = "\t".join(chrinf)
        nullind = ind.replace("NULL",".") 

        varsum_final.write(nullind + "\t.\t.\t.\t.\t.\t.\t.\t.\t"+ hgmd[cleanup] + "\n")








