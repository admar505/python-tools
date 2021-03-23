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
parser.add_argument("--conflicts",help="download from NCBI/clinvar, ",required=True)
#parser.add_argument("--ABthreshold",help="this is a value at which to trust allele ballance, default is 15",default=0.15)
#parser.add_argument("--lc",help="If no genotype can be found, then this is the list of novel or low coverage pages")

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
conflictfi = csv.DictReader(open(args.conflicts,'r'),delimiter='\t')#make a DICT FOR THESE. really much better

varsum_out_name = "hg19_clinvar_" + args.named + ".txt"
varsum_final = open(varsum_out_name,"w")
bigsum_out_name = args.named + "_variant_summary_GRCh37_annovar.txt"
bigsum_final = open(bigsum_out_name,"w")


#---------subroutines-----------------#

def getCLN(varline):#there is a double key, one column has to be split.
                        #this is almost the sole reason there is a key.
    
    varlst = []
    varlst = varline['PhenotypeIDS'].split(",")
    CLNDSDB = {}
    CLNDSDBID = {}


    
    #varlst = varbreak.split("|")#It is unclear as to why the newer clinbar uses commas and 
                                #why it also uses the the "|" I will attempt to use both.

    for varbreak in varlst:
        for clnln in varbreak.split('|'):
      
            array = clnln.split(":")#because clinvar SUCKS, the ID is only the last entry in this record.
            arrct  = 0

            while arrct < (len(array) - 1):
        
                CLNDSDB[array[arrct]] = arrct
                arrct = arrct + 1


            CLNDSDBID[array[-1]] = array[-1]
            ctmp = ",".join(CLNDSDB.keys()) 
            clndsdbid = ctmp.translate(str.maketrans({'|':';'}))


    returnstring = clndsdbid + "\t" + ",".join(CLNDSDBID.keys())
    
    return(returnstring)

def checkGuidelines(glines):#simply checks if guidelines are dashes, it needs to be dots.
    
    if "-" in str(glines):
        glines = "NA"

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
    

def resolveConflicts(diamond,idvar):#$$hashmap of ids to variation, $$current ID
    sig_string = None#initialize as a None so. None.
    load = False 
    bucket = {}##store for uniquifieing and ordering


    if idvar in diamond:##calling it diamonds as it is NOT CONFLICT FREE
        load = True
        for sig1 in diamond[idvar]:
            bucket[sig1] = sig1
            bucket[diamond[idvar][sig1]] = diamond[idvar][sig1] ##note, that diamond[idvar][sig1] SHOULD be SIG2 
   

    if load is True:
        sig_string = " ".join(sorted(bucket.keys()))

    return(sig_string)


#---------------main------------------#

conflicts = {}##NCBI_Variation_ID = {Submitter1_LastEval} = Submitter2_LastEval

for varid in conflictfi:

    if str(varid['NCBI_Variation_ID']) in conflicts:
        conflicts[varid['NCBI_Variation_ID']][varid['Submitter1_ClinSig']] = varid['Submitter2_ClinSig']
    
    else:
        conflicts[varid['NCBI_Variation_ID']] = {} 
        conflicts[varid['NCBI_Variation_ID']][varid['Submitter1_ClinSig']] = varid['Submitter2_ClinSig']


for hln in hgmdfi:#HGMD is going in.

    key = fixStrand(hln)
    hgmd[key] = hln['hgmd_accession'] + "\t"  + hln['variant_class'] + "\t"  + hln['primary_pubmed']


#JUST A HEADER, dont freakout.
varsum_final.write("#Chr\tStart\tEnd\tRef\tAlt\tCLNSIG\tCLNDBN\tCLNACC\tCLNDSDB\tCLNDSDBID\tReview\tGuidelines\tVariationID\thgmdID\thgmdVC\thgmdPUBMED\n")
bigsum_final.write("#Chr\tStart\tEnd\tRef\tAlt\tCLNSIG\tCLNDBN\tCLNACC\tCLNDSDB\tCLNDSDBID\tReview\tGuidelines\trsID\tVariationID\n")





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

    #backtranslate for formats needed.
    changeguidelinetodot = checkGuidelines(var['Guidelines'])
    
    RCVaccession = var['RCVaccession'].translate(str.maketrans({'|':';'}))
    clnDBN = var['RCVaccession'].translate(str.maketrans({'|':';'}))
    clinsig = resolveConflicts(conflicts,var['VariationID'])


    if clinsig is None:
        sigs = var['ClinicalSignificance'] 
        clinsig = sigs.replace(",","") 
         

    varsum_final.write(var['Chromosome'] +"\t"+ var['Start'] + "\t" + var['Stop'] + "\t" + var['ReferenceAlleleVCF'] + "\t" + var['AlternateAlleleVCF']\
             + "\t" + clinsig + "\t" + clnDBN + "\t" + RCVaccession + "\t" + clnDBS + "\t" + var['ReviewStatus']\
             + "\t" + changeguidelinetodot + "\t" + var['#AlleleID'] + "\t" + hgmdinfo  +  "\n")

    bigsum_final.write(var['Chromosome'] +"\t"+ var['Start'] + "\t" + var['Stop'] + "\t" + var['ReferenceAlleleVCF'] + "\t" + var['AlternateAlleleVCF']\
             + "\t" + clinsig + "\t" + clnDBN + "\t" + RCVaccession + "\t" + clnDBS + "\t" + var['ReviewStatus']\
             + "\t" + changeguidelinetodot + "\t" + var['RS# (dbSNP)'] + "\t" +  var['#AlleleID']  +  "\n")

##Now, fill in the items which did not have a match in the clinvar files.

for cleanup in hgmd:

    if cleanup not in removehgmd:
        


        chrinf = cleanup.split("-")
        ind = "\t".join(chrinf)
        nullind = ind.replace("NULL",".") 

        varsum_final.write(nullind + "\t.\t.\t.\t.\t.\t.\t.\t.\t"+ hgmd[cleanup] + "\n")








