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

#----------file-handling-defs----------#

def newName(v,f):

    m = re.compile('([a-zA-Z0-9_-]+)\..*')
    try:
        typeofvcf = m.match(v)
        sampleid = m.match(str(os.path.basename(f)))
        filename = sampleid.group(1) + "." + typeofvcf.group(1) + ".RESULTS.txt"
        return filename

    except AttributeError:

        filename = "TEST.RESULTS.txt"
        return filename

#----------file-handling-main---------#

args = parser.parse_args()
answerfi = args.answer #-->  to dict, use CSV
vcffi = args.vcf       #--> standard pyvcf
fullfi = args.fullvcf  #--> to standard pyvcf
combo = args.combo    # this is array, as this can be several
newres = vgr.Writer(open(newName(vcffi,fullfi),"w"))#temp name.
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


def gtCallOfficial(genocall):#return variant-caller asserted GT
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
    print str(ab) + "\tIN " + str(alt) + " CHECK\t GL TO AB  " + str(gl)

    isinWT = False
    isValidAB = False

    for nucl in gl.split(","):
        if nucl == ref:
            isinWT = True

    variant = 0

    while variant < len(ab):
        print str(ab[variant]) + " starting checks GL2AB"

        if float(ab[variant]) > float(args.ABthreshold) or float(ab[variant]) == 0:# CHECK if AB is better than .15 or 
                                                                                   #set to zero for HOMOVARS ()
            print str(alt[variant]) +  "\tVAR inside the AB check\t" + str(gl)     #this may fail if minor allele get ready for that.

            if str(alt[variant]) in str(alt):
                isValidAB = True

        variant += 1

    if (isinWT is True and isValidAB is True) or (isValidAB is True and callao.amIHOMO is True):
        itisgood = True

    return itisgood



def checkQUAL(correctGT,call,var):
    qual = True

    return qual

def raiseFAIL(intendedcall,reason):##proto:callao,"REASON, in text"
    print str(intendedcall.defGT_Dict) + "\t" + str(reason)

def printHetOrHomo():#for printing will direct to homo or het.
    print "het or homo"


def printVAR(callAO,var_fb,answer,truealt):
    #validalt = callAO.
    
    #print truealt + "\t" + callAO.WT + "\t" + str(callAO.defGT_Dict)

    newrecord = vgr.model._Record(var_fb.CHROM,var_fb.POS,var_fb.REF,truealt,{})

    newrecord.INFO['FBGenoType'] = var_fb.REF  + truealt
    newrecord.INFO['FBRefAlleleCount'] = var_fb.INFO['QR'] #might need to split
    newrecord.INFO['VAPOR_URL'] = answer['heturl']
    newrecord.INFO['EFF_HGVS'] = answer['hethgvs']
    


    newres.write_record(newrecord)


def getGoodALT(calldict):#give a dictionary, and it will return only the TRUE
                         #might be needed as there could be two
    validalt = None

    for arr in calldict:
        if calldict[arr] is True:
            validalt = arr

    return str(validalt)

def assignFinalGT(callAO,var_fb,answer):#ok, so, some logic, if the gt gl all work,
    AB = var_fb.INFO['AB']               #and the qqual is above threshold, and the AB is good, give it a blam.
    qual = var_fb.QUAL                   #see if all agree. if so, see if the call is homo, het or wt
                                         ##then, can pull the lines here from the answer bed: homohgvs  homourl hethgvs heturl  wthgvs  wturl
                                         #Also: ensure, with the rsid, that the call is valid.
    truealt = getGoodALT(callAO.defGT_Dict)
    asserted_gt = gtCallOfficial(callAO)


    best_gtGL = getBestGL(callAO.assGT_GL)#this value stores what should be returned. test all against this value.
    print best_gtGL + "\tB"

    if checkGLtoAB(best_gtGL,AB,qual,var_fb.ALT,var_fb.REF,callAO) is True:#DOES all information match the asserted type,

        if checkQUAL(best_gtGL,callAO,var_fb) is True:
            print "qual is good, send to printer"
            print truealt
            printVAR(callAO,var_fb,answer,truealt)

    else:
        truealt = None#reassign what truealt is
        print best_gtGL
        raiseFAIL(callAO,"UNMATCHED_GENOTYPE")#RIGHT NOW RAISE FAIL --> later, assign correct type.


def returnWT(wtcall,variant,answer):##if I have a obvious WT allele, this just kicks it. eventually turn to printer
    #print answer['wthgvs'] + "\t" + answer['wturl']

    wtrecord = vgr.model._Record(variant.CHROM,variant.POS,variant.REF,variant.REF,{})

    wtrecord.INFO['FBGenoType'] = variant.REF + variant.REF
    wtrecord.INFO['FBRefAlleleCount'] = variant.INFO['QR']
    wtrecord.INFO['VAPOR_URL'] = answer['heturl']

    newres.write_record(wtrecord)



def determineCall(varobj,targ): #This will be the beginning of determining the call.
                                #step TWO
                                #---------------------------------------------------
    for variant in varobj:      #should this differentiate between dels and snps? lets see here.
        #print variant.POS       #get call -> assign to this type --> success.
        #print targ["start"]     #if WT+, cut to chase?

        callobj = loadaltdats.detGenoType(variant)

        if callobj.amIWT is  True:
            returnWT(callobj,variant,targ)#just call it done and returnWT+():

        else:

            try:
                #print "CALL IS GOOD " + str(callobj.defGT_Dict)
                #print "here is this " + str(callobj.assGT_GL)
                assignedGT = assignFinalGT(callobj,variant,targ)
                #print assignedGT

            except AttributeError:#eventually do variant failure return to file
                print  "FAILED to get variant for rsid: " + str(targ['rsid'])


#####----------------MAIN--------------####      #####----------------MAIN--------------####

for bed in bedfi:#as csvDictReader

    #FIRST: collapse ithe rsids into possibles, example A-C or A-T for the same.
    #try:#need to pass the specific bed line that is target

        variant = resvcf.fetch(str(bed['chr']),int(bed['start']),int(bed['stop']))#pull variant
        call = determineCall(variant,bed)                                         #send to make call

    #except ValueError:#initiate error checks. here.
     #   print "WARNING:No variant for answerbed regioni " + answerfi + " " + str(bed['rsid'])













