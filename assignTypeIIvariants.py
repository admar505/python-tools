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
parser.add_argument("--lc",help="If no genotype can be found, then this is the list of novel or low coverage pages",default=0.15)

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
lcfi = args.lc

try:
    resvcf = vcf.Reader(open(vcffi,'r'))
    bedfi  = csv.DictReader(open(answerfi,'r'),delimiter='\t')
    recovery = csv.DictReader(open(lcfi,'r'),delimiter='\t')

except (TypeError,NameError) as e:
    print "\n\n\tUSE -h thanks.\n\n"

results = {}#stores the results lines;
#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def returnWT(wtcall,variant,answer):##if I have a obvious WT allele, this just kicks it. eventually turn to printer
    #print answer['wthgvs'] + "\t" + answer['wturl']

    wtrecord = vgr.model._Record(variant.CHROM,variant.POS,variant.REF,variant.REF,{})

    wtrecord.INFO['FBGenoType'] = variant.REF + variant.REF
    wtrecord.INFO['FBRefAlleleCount'] = variant.INFO['RO']
    wtrecord.INFO['EFF_PROT'] = 'NULL_PROT'
    wtrecord.INFO['VAPOR_URL'] = answer['wturl']
    wtrecord.INFO['EFF_HGVS'] = answer['wthgvs']
    wtrecord.INFO['RSID'] = answer['rsid']
    wtrecord.INFO['FBTotalDepth'] = variant.INFO['DP']
    wtrecord.INFO['QUAL'] = wtcall.retQUAL
    wtrecord.INFO['FBReferenceAlleleQ'] = variant.INFO['QR']

    newres.write_record(wtrecord)

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
    returnable_string =  ",".join(sortalleles)

    return returnable_string

def cigarCheck(call):#for now, return just TRUE or FALSE
    reliable = True

    print "IN CIGAR"


    print str(call.full.POS) + "\t IN CIGAR \t" +  str(call.full.INFO['CIGAR'])

    if re.match('.*?M.*?',str(call.full.INFO['CIGAR'])):
        reliable = False

    print reliable

    return reliable


def checkGLtoAB(gl,ab,ql,alt,ref,callao):  #return True if GL and AB vals agree on call.
                                           #approach: verify that GTbyGL is either WT, or, is high AB alt
    itisgood = False                       #Just discovered that sometimes, when vars are colocated, true GL includes the prev.
    isinWT = False                         #vars positions, and not the correct GT. PLAN: adjust here to see,
    isValidAB = False                      #if CIGAR doesnt match the call, dont use the GL.

    is_cigar_valid = cigarCheck(callao)

    for nucl in gl.split(","):

        if (nucl == ref):
            isinWT = True


    variant = 0

    while variant < len(ab):

        if float(ab[variant]) > float(args.ABthreshold) or float(ab[variant]) == 0:# CHECK if AB is better than .15 or
                                                                                   #set to zero for HOMOVARS ()
            if str(alt[variant]) in str(alt):
                isValidAB = True

        variant += 1

    if ((isinWT is True and isValidAB is True) or (isValidAB is True and callao.amIHOMO is True)) and is_cigar_valid is True:
        itisgood = True

    return itisgood



def checkQUAL(correctGT,call,var):
    qual = True

    return qual

def raiseFAIL(intendedcall,answer,reason):##proto:callao,"REASON, in text"
    print str(intendedcall.defGT_Dict) + "\t" + str(reason) + "\t" + str(answer)

def getABCall(callao):
    homo = True

    fb_gt = gtCallOfficial(callao).split(',')

    if fb_gt[0] != fb_gt[1]:
        homo = False

    print str(fb_gt) + "   Official GT call"

    return homo


def printHetOrHomo(callao,ans,alt,trust_gl):#for printing will direct to homo or het.
    ret = []
    hetgt = [alt,callao.WT]

    print alt + " " + str(callao.full.POS) +" " + str(callao.amIHOMO)

    homo = True


    def __big_join__(a,b):
        if len(a) > 1 or len(b) > 1:
            genotype = "/".join(sorted([a,b]))
            print genotype
        else:
            genotype = "".join(sorted([a,b]))
            print genotype

        return genotype


    if trust_gl is False:
        homo = getABCall(callao)


    if homo is True:
        ret.append(ans['homourl'])
        ret.append(ans['homohgvs'])
        ret.append(__big_join__(alt,alt))

    else:
        ret.append(ans['heturl'])
        ret.append(ans['hethgvs'])
        ret.append(__big_join__(alt,callao.WT))

    return ret


def printVAR(callAO,var_fb,answer,truealt,trust_gl):

    newrecord = vgr.model._Record(var_fb.CHROM,var_fb.POS,var_fb.REF,truealt,{})

    newrecord.INFO['FBGenoType'] = printHetOrHomo(callAO,answer,truealt,trust_gl)[2]
    newrecord.INFO['FBRefAlleleCount'] = var_fb.INFO['RO'] #might need to split
    newrecord.INFO['VAPOR_URL'] = printHetOrHomo(callAO,answer,truealt,trust_gl)[0]
    newrecord.INFO['EFF_HGVS'] = printHetOrHomo(callAO,answer,truealt,trust_gl)[1]
    newrecord.INFO['EFF_PROT'] = 'NULL_PROT'
    newrecord.INFO['RSID'] = answer['rsid']
    newrecord.INFO['FBTotalDepth'] = var_fb.INFO['DP']
    newrecord.INFO['QUAL'] = callAO.retQUAL
    newrecord.INFO['FBReferenceAlleleQ'] = var_fb.INFO['QR']
    newres.write_record(newrecord)


def getGoodALT(calldict):#give a dictionary, and it will return only the TRUE
                         #adjustment might be needed as there could be two
    validalt = None

    for arr in calldict:
        if calldict[arr] is True:
            validalt = arr

    return str(validalt)

def sumAB(vcfvar):#checking to see if the AB is noisy
    totalab = 0

    if vcfvar.INFO['AB']:
        for ab in vcfvar.INFO['AB']:
            totalab += float(ab)

    print str(totalab) + "   IN TOTAL AB"

    return totalab

def assignFinalGT(callAO,var_fb,answer):#ok, so, some logic, if the gt gl all work,
    AB = var_fb.INFO['AB']               #and the qqual is above threshold, and the AB is good, give it a blam.
    qual = var_fb.QUAL                   #see if all agree. if so, see if the call is homo, het or wt
                                         ##then, can pull the lines here from the answer bed: homohgvs  homourl hethgvs heturl  wthgvs  wturl
                                         #Also: ensure, with the rsid, that the call is valid.
    truealt = getGoodALT(callAO.defGT_Dict)
    asserted_gt = gtCallOfficial(callAO)

    #print "POTENTIAL_CALL  " + asserted_gt + "  " + truealt + " +  +"  + str(answer) + " " + str(var_fb.POS) + " " + var_fb.REF

    sumofballance = sumAB(var_fb)

    print sumofballance

    best_gtGL = getBestGL(callAO.assGT_GL)#this value stores what should be returned. test all against this value.

    if checkGLtoAB(best_gtGL,AB,qual,var_fb.ALT,var_fb.REF,callAO) is True:#DOES all information match the asserted type,

        if checkQUAL(best_gtGL,callAO,var_fb) is True:
            printVAR(callAO,var_fb,answer,truealt,True)

    elif sumAB(var_fb)  < float(args.ABthreshold) + float(.025):#IF WT WITH NOISE, this will capture,and send to WT printer.

        returnWT(callAO,var_fb,answer)

    elif (sumAB(var_fb) >= float(args.ABthreshold) + float(.025)) or (float(sumAB(var_fb)) == float(0.0)): #Here, if I turn off the GL trust due to the merged CIGAR issue.
        print "through the catcher"
        if checkQUAL(best_gtGL,callAO,var_fb) is True:
            printVAR(callAO,var_fb,answer,truealt,False)

    else:
        truealt = None                        #reassign what truealt is
        raiseFAIL(callAO,answer,"UNMATCHED_GENOTYPE")#RIGHT NOW RAISE FAIL --> later, assign correct type.

def determineCall(varobj,targ): #This will be the beginning of determining the call.
                                #step TWO
                                #---------------------------------------------------
    for variant in varobj:      #should this differentiate between dels and snps? lets see here.
        #print variant.POS       #get call -> assign to this type --> success.
        #print targ["start"]     #if WT+, cut to chase?
        if int(variant.POS) == int(targ['start']) + 1:#protects from off target calls
        #if int(variant.POS) > int(targ['start']) and int(variant.POS) < int(targ['stop']):#protects from off target calls
            callobj = loadaltdats.detGenoType(variant)


            if callobj.amIWT is  True:
                returnWT(callobj,variant,targ)#just call it done and returnWT+():

            else:

                try:
                    #print "CALL IS GOOD " + str(callobj.defGT_Dict) + "\t" +  str(variant.POS)
                    #print "here is this " + str(callobj.assGT_GL) + "\t" +  str(variant.POS)
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













