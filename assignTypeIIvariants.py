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
from collections import defaultdict
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
#header is:chr start   stop    rsid    homohgvs    homourl hethgvs heturl  wthgvs  wturl
parser = argparse.ArgumentParser(description="Assigns to each Type II variant, a vapor URL. if missing, returns error, but keeps going. Also assigns combos. Uses the new RESULTS.txt library!!")
parser.add_argument("--answer",help="the mapping of variant to VAPor page")
parser.add_argument("--vcf",help="VCF file, bgzipped, indexed with tabix -p vcf")
parser.add_argument("--fullvcf",help="FULL genome VCF file, bgzipped, indexed with tabix -p vcf")
parser.add_argument("--combo",help="combination types",action='append')
parser.add_argument("--ABthreshold",help="this is a value at which to trust allele ballance, default is 15",default=0.15)
parser.add_argument("--lc",help="If no genotype can be found, then this is the list of novel or low coverage pages",default=0.15)
parser.add_argument("--hap",help="vcf-allele-haplotypes output",required=True)
parser.add_argument("--trans",help="the gene to transcript file for PGX")
parser.add_argument("--rpt",help='the repeats formatted file, reuse argument for multiple searches.\n format is: \"chrom start stop base leader unit wt_version wt_length hgvs_primitive\"',action='append')
parser.add_argument("--add",help="inject a line, usually used to add results from an external tool and add it to the combo search sections",action='append')

all_hgvs = {}#this will store all the hgvs for combo type.

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
hapfi = args.hap
addfiles = args.add
pgx_transfi = args.trans
repeats = args.rpt

try:
    resvcf = vcf.Reader(open(vcffi,'r'))
    fullvcf = vcf.Reader(open(fullfi,'r'))
    bedfi  = csv.DictReader(open(answerfi,'r'),delimiter='\t')
    recovery = open(lcfi,'r')
    haps = csv.DictReader(open(hapfi,'r'),delimiter=',')
    pgx_trans = open(pgx_transfi,'r')

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
    all_hgvs[str(answer['wthgvs'])] = str(answer['rsid'])
    all_hgvs[str(variant.REF + variant.REF)] = str(answer['rsid'])

    newres.write_record(wtrecord)

def getBestGL(gtgl):#return best GT given GL
    best = None

    sortedgtgl = sorted(gtgl, key=gtgl.get, reverse=True)
    best = sortedgtgl[0]

    return best


def gtCallOfficial(genocall):#return variant-caller asserted GT, for when GL is unreliable
    alleles = []             #

    sortalleles = sorted(genocall.get_default_GT)
    returnable_string =  ",".join(sortalleles)

    #print str(alleles) + "\tIN GT CALL OFFICIAL\t"  + str(genocall.full.POS)

    return returnable_string

def cigarCheck(call):#for now, return just TRUE or FALSE
    reliable = True

    if re.match('.*?M.*?',str(call.full.INFO['CIGAR'])):
        reliable = False

    #print str(reliable) + "\t" +  str(call.full.POS) + "\t IN CIGAR \t" +  str(call.full.INFO['CIGAR'])

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


def checkQUAL(correctGT,call,var):#this assumes that the highest
    qual = False

    if var.QUAL >= 100:
        qual = True

    if int(sorted(var.INFO['QR'],reverse=True)[0]) >= 100:
        qual = True

    return qual


def raiseFAIL(intendedcall,recover,reason,answer,types):##proto:callao,"REASON, in text"
                                                        #types is either ":LC" for normal SNPs,
                                                        #or "Unresolved"
    def __get_lc__(rec,line_selection,typer):
        lc_return = None

        recs = csv.DictReader(rec,delimiter='\t')
        rec.seek(0)

        hgvs = line_selection.split(':')

        specific = str(hgvs[1]) + str(typer) + "||" + str(hgvs[0])

        for recline in recs:
            if str(recline['failhgvs']) == str(specific):
                lc_return = str(recline['failurl'])
        return lc_return


    failrecord = vgr.model._Record(intendedcall.full.CHROM,intendedcall.full.POS,intendedcall.full.REF,intendedcall.full.REF,{})

    failrecord.INFO['FBGenoType'] = "Unresolved"
    failrecord.INFO['FBRefAlleleCount'] = intendedcall.full.INFO['RO']
    failrecord.INFO['EFF_PROT'] = 'NULL_PROT'
    failrecord.INFO['EFF_HGVS'] = str(answer['homohgvs']) + str(types)
    failrecord.INFO['VAPOR_URL'] = __get_lc__(recover,answer['homohgvs'],types)
    failrecord.INFO['RSID'] = answer['rsid']#dont forget to trim this biotch
    failrecord.INFO['FBTotalDepth'] = intendedcall.full.INFO['DP']
    failrecord.INFO['QUAL'] = intendedcall.full.QUAL
    failrecord.INFO['FBReferenceAlleleQ'] = intendedcall.full.INFO['QR']
    failrecord.INFO['FAIL_REASON'] = reason

    all_hgvs[str("Unresolved")] = answer['rsid']
    all_hgvs[str(answer['homohgvs']) + str(types)] = answer['rsid']

    newres.write_record(failrecord)


def getABCall(callao):
    homo = True

    fb_gt = gtCallOfficial(callao).split(',')

    if fb_gt[0] != fb_gt[1]:
        homo = False

    #print str(fb_gt) + "   Official GT call" + str(callao.full.POS)

    return homo


def check_final_GT(callao,trustgl):#return GL call if trust is true, else use the AB

    genotype = None
    if trustgl is True:
        #print str(callao.get_default_GT) + "\tTRUST IS TRUE\t" + str(callao.defGT_Dict)
        genotype = callao.get_default_GT

    else:
        #print str(callao.defGT_Dict) + "\tTRUST is FALSE\t" + str(callao.get_default_GT)

        genotype = []
        for val in callao.defGT_Dict:

            #print callao.defGT_Dict[val]
            if callao.defGT_Dict[val] is True:
                genotype.append(str(val))

    #print genotype
    return sorted(genotype)

def check_final_Homo(gt):

    homo = True

    if gt[0] != gt[1]:
        homo = False

    return homo

def check_final_Het(genotype, answered):
    het = False
    if str(''.join(sorted(genotype))) == str(answered):
        het = True

    return het


def choose_answer(homo,alt,callao,answer,trust_gl):#THIS is final check, needs to make sure, that if
                                          #this is het, then the het is ref/alt, not alt1/alt2,
                                          #and if homo, not homo alt2/alt2 that is not covered by the
                                          #answer bed.

    the_real_gt = check_final_GT(callao,trust_gl)

    allowed_answer_line = None

    for ans_alt in  answer['calls']:#CALLS contains the ALT.
        #print str(the_real_gt) + "\tin choose answer, the check for the real GT"
        if homo is True and check_final_Homo(the_real_gt) is True:  #if HOM, then both the alt matches alt
                                                                    #and the realGT is Homo

            if alt == ans_alt:#here should be checkFinalGT

                allowed_answer_line = answer['calls'][ans_alt]
                #print "captured HOMO " + str(answer['calls'][ans_alt]['homourl'])

        else:

            for potential_call in  callao.full.ALT:
                if potential_call == ans_alt and check_final_Het(the_real_gt,''.join(sorted(callao.full.REF,potential_call))) is True:#here should be checkFinalGT

                    allowed_answer_line = answer['calls'][ans_alt]
                    #print "captured HET  " + str(potential_call) + "\t" +  str( answer['calls'][ans_alt]['heturl'])

    if allowed_answer_line is None:
        raiseFAIL(callao,recovery,"UNKNOWN_GENOTYPE",answer['wt'],":LC")
        #print "CAPTURED FAIL    "  + str(answer)

    else:
        return allowed_answer_line


def printHetOrHomo(callao,ans,alt,trust_gl):#for printing will direct to homo or het.
    ret = []                                #goal, decide if gl is trustable, if so, use it, if not
                                            #use AB

    #print alt + " " + str(callao.full.POS) +" in het or   "  +  str(trust_gl)   +   " homo decider " + str(callao.amIHOMO)

    homo = True

    def __big_join__(a,b):#to format the gt correctly, with slash or not.
        if len(a) > 1 or len(b) > 1:
            genotype = "/".join(sorted([a,b]))
        else:
            genotype = "".join(sorted([a,b]))

        return genotype


    if trust_gl is False:
        homo = getABCall(callao)
        #print str(homo) + " " + str(callao.full.POS)
    elif callao.amIHOMO is False:
        homo = False

    #I think, here I can check, and see if it can be redirected to reportFAIL. test with also the undet type.
    #
    final_ans = choose_answer(homo,alt,callao,ans,trust_gl)
    #print str(final_ans) + "\tWHAT IS THE FINAL ANSWER"

    if homo is True and final_ans is not None:#final_ans could return None if no good answer is available.
        ret.append(final_ans['homourl'])
        ret.append(final_ans['homohgvs'])
        ret.append(__big_join__(alt,alt))

    elif final_ans is not None:
        ret.append(final_ans['heturl'])
        ret.append(final_ans['hethgvs'])
        ret.append(__big_join__(alt,callao.WT))


    return ret

def printVAR(callAO,var_fb,answer,truealt,trust_gl):
    if len(printHetOrHomo(callAO,answer,truealt,trust_gl)) > 0:#protects against no good answer

        newrecord = vgr.model._Record(var_fb.CHROM,var_fb.POS,var_fb.REF,truealt,{})

        newrecord.INFO['FBGenoType'] = printHetOrHomo(callAO,answer,truealt,trust_gl)[2]
        newrecord.INFO['FBRefAlleleCount'] = var_fb.INFO['RO'] #might need to split
        newrecord.INFO['VAPOR_URL'] = printHetOrHomo(callAO,answer,truealt,trust_gl)[0]
        newrecord.INFO['EFF_HGVS'] = printHetOrHomo(callAO,answer,truealt,trust_gl)[1]
        newrecord.INFO['EFF_PROT'] = 'NULL_PROT'
        newrecord.INFO['RSID'] = answer['wt']['rsid']
        newrecord.INFO['FBTotalDepth'] = var_fb.INFO['DP']
        newrecord.INFO['QUAL'] = callAO.retQUAL
        newrecord.INFO['FBReferenceAlleleQ'] = var_fb.INFO['QR']

        all_hgvs[printHetOrHomo(callAO,answer,truealt,trust_gl)[2]] =  answer['wt']['rsid']
        all_hgvs[printHetOrHomo(callAO,answer,truealt,trust_gl)[1]] =  answer['wt']['rsid']

        newres.write_record(newrecord)

    #eventually put in fail here. doing that now though? it does kick the fail. and just doesnt print here.

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

    #print str(totalab) + "   IN TOTAL AB"

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
    best_gtGL = getBestGL(callAO.assGT_GL)#this value stores what should be returned. test all against this value.

    if checkGLtoAB(best_gtGL,AB,qual,var_fb.ALT,var_fb.REF,callAO) is True:#DOES all information match the asserted type,

        #print str(sumofballance) + "\tTOP OF SELECTION\t" +  str(var_fb.POS)
        if checkQUAL(best_gtGL,callAO,var_fb) is True:
            printVAR(callAO,var_fb,answer,truealt,True)
        else:
            raiseFAIL(callAO,recovery,"LOW_QUALITY_GENOTYPE",answer['wt'],":LC")

    elif sumAB(var_fb)  < float(args.ABthreshold) + float(.025) and float(sumAB(var_fb)) != float(0.0):#IF WT WITH NOISE, this will capture,and send to WT printer.

       # print str(sumofballance) + "\tNOISE REDIRECT ELIF\t" +  str(var_fb.POS)

        if checkQUAL(best_gtGL,callAO,var_fb) is True:
            returnWT(callAO,var_fb,answer['wt'])

        else:
            raiseFAIL(callAO,recovery,"LOW_QUALITY_GENOTYPE",answer['wt'],":LC")

    elif (sumAB(var_fb) >= float(args.ABthreshold) + float(.025)) or (float(sumAB(var_fb)) == float(0.0)): #Here, if I turn off the GL trust due to the merged CIGAR issue.
       # print "through the catcher\t" + str(sumofballance) + "\t" +  str(var_fb.POS) + "\t" + truealt
        if checkQUAL(best_gtGL,callAO,var_fb) is True:

            printVAR(callAO,var_fb,answer,truealt,False)
        else:
            raiseFAIL(callAO,recovery,"LOW_QUALITY_GENOTYPE",answer['wt'],":LC")

    else:
        truealt = None                        #reassign what truealt is
        raiseFAIL(callAO,recovery,"UNMATCHED_GENOTYPE",answer['wt'],":LC")#RIGHT NOW RAISE FAIL --> later, assign correct type.


def determineCall(varobj,rsindex,all_bed): #This will be the beginning of determining the call.
                                #step TWO
                                #---------------------------------------------------
    for variant in varobj:      #should this differentiate between dels and snps? lets see here.
        #print variant.POS       #get call -> assign to this type --> success.
        #print targ["start"]     #if WT+, cut to chase?


        try:

            if int(variant.POS) == int(all_bed[rsindex]['start']) + 1:#protects from off target calls
        #if int(variant.POS) > int(targ['start']) and int(variant.POS) < int(targ['stop']):#protects from off target calls
               #targ = choose_answer(all_bed,variant,rsindex)
                callobj = loadaltdats.detGenoType(variant)

                if callobj.amIWT is  True:
                    returnWT(callobj,variant,all_bed[rsindex]['wt'])#just call it done and returnWT+():

                else:

                    #print "CALL IS GOOD " + str(callobj.defGT_Dict) + "\t" +  str(variant.POS)
                    #print "here is this " + str(callobj.assGT_GL) + "\t" +  str(variant.POS)
                    assignedGT = assignFinalGT(callobj,variant,all_bed[rsindex])
                    #print assignedGT

        except AttributeError:#eventually do variant failure return to file

                print  "FAILED to get variant for rsid: " + str(rsindex)


####---------------PRE_LOAD_DEFS----####

def tree(): return defaultdict(tree)


def answer_dict(full_ans):#create tree of rsid to answer.

    ret_dict = tree()
    #ret_dict = defaultdict(tree)

    #$lets go different, like, a sub that detects multiple, in one pass, then, goes through and decide what in second?
    for line in full_ans:
        rsnum = re.search('(rs\d+)\-(\w+)\-(\w+)',str(line['rsid']))

        ret_dict[rsnum.group(1)]['calls'][rsnum.group(3)] = line
        ret_dict[rsnum.group(1)]['start'] = line['start']
        ret_dict[rsnum.group(1)]['stop'] = line['stop']
        ret_dict[rsnum.group(1)]['chr'] = line['chr']
        ret_dict[rsnum.group(1)]['wt'] = line

    return ret_dict

#############-----------------HAPloTYPE_SUb_and-DEFS----------------------###########
#overall goal:for gene, send each element to hapformatter, return type. save type in gene specific


def hapFormat(haptype_inc):#BLAH, I think GEne specific regex. darn.well, break em out by bite, the send to recon
                       #method, and return the

    haptype = []


    try:#
        hid = re.search('(\w+)(\*\w)(\w?)(\w?)(\(?[A-Za-z,]{0,10}\)?)',haptype_inc)

        if hid.group(1) == "TPMT":      #TPMT, detect if hid contains subtype, and allow
            haptype.append(hid.group(2))#

            if hid.group(3).isalpha() is True and hid.group(3).isdigit() is True:
                haptype.append(hid.group(3))

            if hid.group(4).isalpha() is True and hid.group(4).isdigit() is True:
                haptype.append(hid.group(4))

        else:
            haptype.append(hid.group(2))

            if hid.group(3).isdigit() is True:
                haptype.append(hid.group(3))

            if hid.group(4).isdigit() is True:
                haptype.append(hid.group(4))

    except TypeError:#
        print "unable to parse " + str(haptype) + " correctly, please check output and code"

    return "".join(haptype)


def mapHaps(gene_name,hapdct):#return collapsed formatted diplotype, two way also.

    formatted = {}

    for haptype in hapdct[gene_name]:
        #print haptype
                                         #check for None,
        formatted[hapFormat(haptype)] = 1#load the haplotype, formatted as a key.

    return formatted


def searchHaps(hapa,hapb,rec,pgxs_handle,genesym):
    #print hapa
    hgvs_and_url = {}

    def __loadtrans__(pfi):
        ret_dict = {}

        for pf in pfi:
            ret_dict[pf['sym']] = pf['trans']

        return ret_dict


    recs = csv.DictReader(rec,delimiter='\t')#I do it like this because I have to rewind and reuse alot. fuck.py for making em like this.
    rec.seek(0)

    pgxs = __loadtrans__(csv.DictReader(pgxs_handle,delimiter='\t'))
    pgxs_handle.seek(0)
    hapmap = {}#

    for line in recs:                               #loading the hapmap, make sure in the lookup table that
        hapmap[line['failhgvs']] = line['failurl']  #all the subtypes that need to be removed ARE.

    hgvslinea = ""
    hgvslineb = ""

    #next, construct the hgvs.. recon and search
    if hapa == "Unresovled":
        hgvslinea  = hapa + "||" + pgxs[genesym]
    else:
        hgvslinea = hapa + "/" + hapb + "||" + pgxs[genesym]
        hgvslineb = hapb + "/" + hapa + "||" + pgxs[genesym]

    #print hgvslinea + "\t" + hgvslineb

    if hgvslinea in hapmap.keys():
        hgvs_and_url['wthgvs'] = hgvslinea
        hgvs_and_url['wturl'] = hapmap[hgvslinea]


    elif hgvslineb in hapmap.keys():
        hgvs_and_url['wthgvs'] = hgvslineb
        hgvs_and_url['hgvslineb'] = hapmap[hgvslineb]

    return hgvs_and_url


def printHap(pgxs_handle,vap_url,genesym):

    def __loadpgx__(pfi,genesym):
        pgxd = {}
        for pf in pfi:
            if pf['sym'] == genesym:
                pgxd = pf

        return pgxd

    def __gtonly__(string):
        gtwood =string.split('|')
        return gtwood[0]


    def __fliphgvs__(wth):

        if re.match('.*?\|\|.*',wth):
            flp = re.search('(.*?)\|\|(.*)',wth)
            wth = flp.group(2) + ":" + flp.group(1)
        return wth

    pgxs = __loadpgx__(csv.DictReader(pgxs_handle,delimiter='\t'),genesym)
    pgxs_handle.seek(0)

    newrecord = vgr.model._Record(pgxs['chr'],pgxs['start'],"NULL","NULL",{})
    newrecord.INFO['FBGenoType'] = __gtonly__(vap_url['wthgvs'])
    newrecord.INFO['VAPOR_URL'] = vap_url['wturl']
    newrecord.INFO['EFF_HGVS'] = __fliphgvs__(vap_url['wthgvs'])
    newrecord.INFO['EFF_PROT'] = 'NULL_PROT'
    newrecord.INFO['RSID'] = 'Haplotype'

    all_hgvs[vapurl['wthgvs']] = 'Haplotype'
    all_hgvs[__gtonly__(vap_url['wthgvs'])] = 'Haplotype'

    newres.write_record(newrecord)

#####----------------REPEAT-DEFS-------####     ######----------------------------------####

def searchRepeat(rpt,rec,s2t):

    sym = rpt.rptSYM
    vapurl = {}

    def __loadtrans__(pfi):
        ret_dict = {}

        for pf in pfi:
            ret_dict[pf['sym']] = pf['trans']

        return ret_dict

    recs = csv.DictReader(rec,delimiter='\t')
    rec.seek(0)

    genetrans = __loadtrans__(csv.DictReader(s2t,delimiter='\t'))
    s2t.seek(0)

    repeat_var = rpt.countRepeats[1] + "||" + str(genetrans[str(sym)])
    vapurl['wthgvs'] = "LC||" + str(genetrans[str(sym)])#I am changing failhgvs to wthgvs to reuse the happrinter.
    #first, load failure contingency

    for recfail in recs:
        if str(recfail['failhgvs']) == vapurl['wthgvs']:
            vapurl['wturl'] = recfail['failurl']

    rec.seek(0)

    #then grab the regular line.
    for recline in recs:

        if str(recline['failhgvs']) == str(repeat_var):

            vapurl['wthgvs'] = str(repeat_var)
            vapurl['wturl'] = recline['failurl']

    return vapurl



def reportRepeat(vcfs,reps,lookuptab,sym2trans):#vcfull,repeat_specs,allthe urls,sym2transand pos file

    repeat = loadaltdats.detRepeats(vcfs,reps)
    #print repeat.countRepeats
    vapurl = searchRepeat(repeat,lookuptab,sym2trans)

    printHap(pgx_trans,vapurl,"UGT1A1")


#####------------INJECT-DEFS-----------###       ######---------INJECT-DEFS-------------###



def addRes(addline):

    def __returnVal__(adder,val):
        for potential in adder:
            if re.match('.*?=.*',potential):
                groups = re.match('(\S+?)\=(\S+)',potential)
                if groups.group(1) == val:
                    return groups.group(2)

    adder = addline.split("\t")
    #print adder
    newrecord = vgr.model._Record(adder[0],1000,"NULL","NULL",{})
    newrecord.INFO['FBGenoType'] =  __returnVal__(adder,'FBGenoType')
    newrecord.INFO['VAPOR_URL'] = __returnVal__(adder,'VAPOR_URL')
    newrecord.INFO['EFF_HGVS'] =__returnVal__(adder,'EFF_HGVS')
    newrecord.INFO['EFF_PROT'] = 'NULL_PROT'
    newrecord.INFO['RSID'] = __returnVal__(adder,'RSID')
    all_hgvs[__returnVal__(adder,'FBGenoType')] = __returnVal__(adder,'RSID')
    all_hgvs[__returnVal__(adder,'EFF_HGVS')] = __returnVal__(adder,'RSID')

    newres.write_record(newrecord)

#####---------combo-value-meals--------####     #####----------------combo-menus--------####

def printCombo(adder,gt,effhgvs,vapurl,rsid):#I feel like a failure in many ways, I never had a way of making the printer generic.
                #I will work on this in the future.

    newrecord = vgr.model._Record(adder,1000,"NULL","NULL",{})
    newrecord.INFO['FBGenoType'] = gt
    newrecord.INFO['VAPOR_URL'] = vapurl.strip()
    newrecord.INFO['EFF_HGVS'] = effhgvs
    newrecord.INFO['EFF_PROT'] = 'NULL_PROT'
    newrecord.INFO['RSID'] = rsid
    newres.write_record(newrecord)

def testCombo(pvars,hgvs):
    test = False

    var_arr = pvars.split(',')
    success = len(var_arr)

    rsids = []

    for variant in var_arr:
        #print variant
        if variant in hgvs.keys():
            success -= 1
            rsids.append(hgvs[variant])

    if success == 0:
        test = True

    return test,','.join(rsids)


def assignCombo(combo,hgvs_list,filename):#approach:break combo down in to the various components, and search for each?
                                 #also, have fail contingency for not found.just a count and success strategy.

    def __gtmaker__(variants):
        gt = []
        var_arr = variants.split(',')
        for gttype in var_arr:
            gtt = gttype.split(':')
            gt.append(gtt[1])

        return ','.join(gt)

    def __filenamer__(fname):
        fn = fname.split('.')
        return str(os.path.basename(fn[0])) + "_COMBO_NOT_FOUND"


    captured = 0
    for possible_combo in combo:
        (variants,url) = possible_combo.split('\t')
        combogood,rsid = testCombo(variants,hgvs_list)

        if combogood is True:
            printCombo('chrN',__gtmaker__(variants),variants,url,rsid)
            captured = 1

    if captured == 0:
        printCombo('chrN',__filenamer__(filename),"NOT_FOUND",url,"None")


#####----------------MAIN--------------####     #####----------------MAIN--------------####

bed_dict= answer_dict(bedfi)

for rsindex in bed_dict:#as csvDictReader

    #FIRST: collapse ithe rsids into possibles, example A-C or A-T for the same.
    #try:#need to pass the specific bed line that is target
    #print rsindex

    variant = resvcf.fetch(str(bed_dict[rsindex]['chr']),int(bed_dict[rsindex]['start']),int(bed_dict[rsindex]['stop']))#pull variant
    #print variant#compatibility test
    call = determineCall(variant,rsindex,bed_dict)                  #send to make call

    #except ValueError:#initiate error checks. here. SEND to checker for
       # print "WARNING:No variant for answerbed regioni " + answerfi + " " + str(bed['rsid'])

hapdat = {} #NOW, get the haplotypes worked out.

for haplotype in haps:#loading patient haps in.
    #print haplotype
    if haplotype['gene'] not in hapdat.keys():
        hapdat[haplotype['gene']] = list()

    hapdat[haplotype['gene']].append(haplotype['allele1'])
    hapdat[haplotype['gene']].append(haplotype['allele2'])

#haplotype handling, contains methods for discerning diplotypes and hap|haplotypes.
for gene_id in hapdat:
    typed = mapHaps(gene_id,hapdat)#I would like to send to printer from 'ere, or fail from 'ere
    hap_keys = typed.keys()
    vapfailurl = searchHaps('Unresovled','',recovery,pgx_trans,gene_id)#if more than two, send to raiseFail

    if len(hap_keys) == 1: #counter, so here, if one, make diploptype homozygous.
        vapurl = searchHaps(hap_keys[0],hap_keys[0],recovery,pgx_trans,gene_id)#if more than two, send to raiseFail

        if len(vapurl.keys()) == 0:#check for empty
            printHap(pgx_trans,vapfailurl,gene_id)

        else:
            printHap(pgx_trans,vapurl,gene_id)

    elif len(hap_keys) == 2:
        vapurl = searchHaps(hap_keys[0],hap_keys[1],recovery,pgx_trans,gene_id)#

        if len(vapurl.keys()) == 0:#check for empty
            printHap(pgx_trans,vapfailurl,gene_id)#HERE, do this:  make this make a normal looking file.

        else:
            printHap(pgx_trans,vapurl,gene_id)

    else:
        printHap(pgx_trans,vapfailurl,gene_id)
#injection, to add a line to be included in the combos (see star5 handling)
if addfiles:

    for addfile in addfiles:
        addres = open(addfile,'r')

        for adds in addres:
            addRes(adds)

#repeats handling:
if repeats:
    for repeat in repeats:
        repeat_set = open(repeat,'r')
        repeat_capture = reportRepeat(fullvcf,repeat_set,recovery,pgx_trans)

#COMBO_types#TRY to reuse the hapPrinter, it would save a lot of code.
if combo:
    for combs in combo:
        combo_set = open(combs,'r')
        assignCombo(combo_set,all_hgvs,combs)



#for hgvs in global_captured_hgvs:







#current full command:assignTypeIIvariants.py --answer PGX.ans.test.bed  --vcf pr.UNK.PGX.good.vcf.gz   --fullvcf 160406_S2_combined.vcf.gz     --combo /ref/NC_PGX/cftrwt.combo.scf --ABthreshold .15 --lc  /ref/NC_PGX/LC.Novel.HaploType.lst  --hap Sample_160406S4-EXT-ERRORs-alleles.csv  --trans /ref/NC_PGX/pgx.trans.lst  --rpt /ref/NC_PGX/UGT1A1_TA.rpt --add combo.tmp




























