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


def assignFinalGT(callfb,AB,qual):#ok, so, some logic, if the gt gl all work,
                                         #and the qqual is above threshold, and the AB is good, give it a blam.
                                         #see if all agree. if so, see if the call is homo, het or wt
                                         ##then, can pull the lines here from the answer bed: homohgvs  homourl hethgvs heturl  wthgvs  wturl
                                         #Also: ensure, with the rsid, that the call is valid.
    asserted_gt = gtCallOfficial(callfb)

    best_gtGL = getBestGL(callfb.assGT_GL)





    for vals in callfb.defGT_Dict:
        print str(callfb.defGT_Dict[vals]) +"\tVALS\t"+ str(vals) + "\tbest\t"  + best_gtGL + "\t" + str(AB) + " assertedGT " +  str(asserted_gt)





def returnWT(wtcall,answerline):
    print answerline['wthgvs'] + "\t" + answerline['wturl']

def determineCall(varobj,targ): #this will be the beginning of determining the call.



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
                assignedGT = assignFinalGT(callobj,variant.INFO['AB'],variant.QUAL)


            except AttributeError:
                print  "FAILED to get variant for rsid: " + str(targ['rsid'])






#####----------------MAIN--------------####      #####----------------MAIN--------------####

for bed in bedfi:#as csvDictReader

    #FIRST: collapse ithe
    #try:#need to pass the specific bed line that is target

        variant = resvcf.fetch(str(bed['chr']),int(bed['start']),int(bed['stop']))
        call = determineCall(variant,bed)

    #except ValueError:#initiate error checks. here.
     #   print "WARNING:No variant for answerbed regioni " + answerfi + " " + str(bed['rsid'])













