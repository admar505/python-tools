#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/home/nucleo/lib/PyVGRes')
sys.path.append('/home/nucleo/lib/AltObject')
import vgr
import altobject
import loadaltdats
import csv
from subprocess import call
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
#chr start   stop    rsid    homohgvs    homourl hethgvs heturl  wthgvs  wturl
parser = argparse.ArgumentParser(description="Assigns to each Type II variant, a vapor URL. if missing, returns error, but keeps going. Also assigns combos. Uses the new RESULTS.txt library!!")
parser.add_argument("--answer",help="the mapping of variant to VAPor page")
parser.add_argument("--vcf",help="VCF file, bgzipped, indexed with tabix -p vcf")
parser.add_argument("--fullvcf",help="FULL genome VCF file, bgzipped, indexed with tabix -p vcf")
parser.add_argument("--combo",help="combination types",action='append')
parser.add_argument("--ABthreshold",help="this is a value at which to trust allele ballance, default is 15",default=0.15)



args = parser.parse_args()
answerfi = args.answer #-->  to dict, use CSV
vcffi = args.vcf       #--> standard pyvcf
fullfi = args.fullvcf  #--> to standard pyvcf
combo = args.combo    # this is array, as this can be several
#special load for vcffi? (vcffi,"r")
newres = vgr.Writer(open('REQUIRED.RESULTS.txt',"w"))#temp name.
try:
    resvcf = vcf.Reader(open(vcffi,'r'))
    bedfi  = csv.DictReader(open(answerfi,'r'),delimiter='\t')
except (TypeError,NameError) as e:
    print "\n\n\tUSE -h thanks.\n\n"

results = {}#stores the results lines;
#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

def getBestGL(gtgl):#return best GT given GL
    best = None

    sortedgtgl = sorted(gtgl, key=gtgl.get, reverse=True)
    best = sortedgtgl[0]

    return best


def gtCallOfficial(genocall):#return variant-caller asserted GT #IAMHERE
    alleles = []

    for alt in genocall.defGT_Dict:

        if genocall.defGT_Dict[alt] is True:
            alleles.append(str(alt))

    if len(alleles) == 1:
        alleles.append(str(genocall.retREF))

    sortalleles = sorted(alleles)
    returnable_string =  "".join(sortalleles)

    return returnable_string


def checkGLtoAB(gl,ab,ql,alt,ref,callao): #return True if GL and AB vals agree on call.
                                   #approach: verify that GTbyGL is either WT, or, is high AB alt
    itisgood = False
    #print str(ab) + "\t" + str(alt) + "\t" + str(gl)

    isinWT = False
    isValidAB = False

    for nucl in gl.split(","):
        if nucl == ref:

            isinWT = True

    variant = 0

    while variant < len(ab):

        if float(ab[variant]) > float(args.ABthreshold):
            print str(alt[variant]) +  "\tis this in the variant and everything here\t" + str(gl)

            if str(alt[variant]) in str(alt):
                isValidAB = True

        variant += 1

    if (isinWT is True and isValidAB is True) or (isValidAB is True and callao.amIHOMO is True):
        itisgood = True

    return itisgood

def assignFinalGT(callAO,var_fb):#ok, so, some logic, if the gt gl all work,
    AB = var_fb.INFO['AB']               #and the qqual is above threshold, and the AB is good, give it a blam.
    qual = var_fb.QUAL                   #see if all agree. if so, see if the call is homo, het or wt
                                         ##then, can pull the lines here from the answer bed: homohgvs  homourl hethgvs heturl  wthgvs  wturl
                                         #Also: ensure, with the rsid, that the call is valid.
    asserted_gt = gtCallOfficial(callAO)

    best_gtGL = getBestGL(callAO.assGT_GL)#this value stores what should be returned. test all against this value.

    if checkGLtoAB(best_gtGL,AB,qual,var_fb.ALT,var_fb.REF,callAO) is True:
        print "it is"
        print str(callAO.amIHOMO) + "\tHOMO:"



    for vals in callAO.defGT_Dict:
        print str(callAO.defGT_Dict[vals]) +"\tVALS\t"+ str(vals) + "\tbest\t"  + best_gtGL + "\t" + str(AB) + " assertedGT " +  str(asserted_gt)



    return best_gtGL

def returnWT(wtcall,answerline):
    print answerline['wthgvs'] + "\t" + answerline['wturl']

def determineCall(varobj,targ): #This will be the beginning of determining the call.
                                #step TWO


    for variant in varobj:      #should this differentiate between dels and snps? lets see here.
        #print variant.POS       #get call -> assign to this type --> success.
        #print targ["start"]     #if WT+, cut to chase?

        callobj = loadaltdats.detGenoType(variant)

        if callobj.amIWT is  True:
            returnWT(callobj,targ)#just call it done and returnWT+():

        else:

            try:
                print "CALL IS GOOD " + str(callobj.defGT_Dict)
                print "here is this " + str(callobj.assGT_GL)
                assignedGT = assignFinalGT(callobj,variant)


            except AttributeError:
                print  "FAILED to get variant for rsid: " + str(targ['rsid'])






#####----------------MAIN--------------####      #####----------------MAIN--------------####

for bed in bedfi:#as csvDictReader

    #FIRST: collapse ithe rsids into possibles, example A-C or A-T for the same.
    #try:#need to pass the specific bed line that is target

        variant = resvcf.fetch(str(bed['chr']),int(bed['start']),int(bed['stop']))
        call = determineCall(variant,bed)

    #except ValueError:#initiate error checks. here.
     #   print "WARNING:No variant for answerbed regioni " + answerfi + " " + str(bed['rsid'])













