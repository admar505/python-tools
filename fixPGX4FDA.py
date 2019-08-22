#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv

parser = argparse.ArgumentParser(description="for TypeII RESULTS files, compare new with old. thanks")
parser.add_argument("--vapor",help="Current VAPOR dump",required=True)
parser.add_argument("--map",help="the file for translating the Metabolizer status",required=True)
parser.add_argument("--answer",help="Answer file, maybe. ")
parser.add_argument("--trans",help="transcript to name mapping")
args = parser.parse_args()
oldfi = args.vapor
newfi = args.map
ansfi = args.answer
trnsfi = args.trans


mapfi = open(newfi,'r')
maps = csv.DictReader(mapfi,delimiter="\t")

ansfi = open(ansfi,'r')
answer = csv.DictReader(ansfi,delimiter="\t")

vaporfi =  open(oldfi,'r')
vapors =  csv.DictReader(vaporfi)
csvfile = oldfi + ".new.csv"
newvapor = csv.writer(open(csvfile, 'w'),  delimiter=',',quoting=csv.QUOTE_ALL)
transf = open(trnsfi,'r')
trans = csv.DictReader(transf,delimiter="\t")

#parse results in a map or dict, or what??
#vapors: is in a dict.
#maps: construct a look up table of of var-> for snps:
####map rsid and var to hgvs
####    match the rsid to rsid, and match the ref and alt, should match, if not found, DIE.
####    add in lookup table.
####    give the proposed UUID is key, value is the Metabolizer status
####readthrough the file and check the Metabolizer Status for each.
####    if good, do nothing.
####    if not good or blank, update file.
####        emit output <--> the update that was made

#-------------------------------------here by DEFSgONS!!----------------------------------*

def getVal(desiredval,answer):#the line being needed to map == desiredval, the full answerfile == answer
    variant_uuid = None

    def __match_rs__(rs_targ):#get correct rsid
        ansfi.seek(0)

        def __split__(nucl):#decide if snp or del
            if len(nucl) == 2:
                return(list(nucl))

            else:
                alts = nucl

                if "/" in nucl:##find the right splitters.
                    alts = nucl.split("/")#here, I am getting Del worked out

                elif "-" in nucl:
                    alts = nucl.split("-")#suggested here

                if 'D' in alts[0]:
                    return(alts[0].lower(),alts[1].lower())

                elif 'D' in alts[1]:
                    return(alts[0],alts[1].lower())

                else:
                    return(alts[0],alts[1])

        def __matchalt__(targ,ansl,nucl,full_answer,desval):#get the right hgvs nucleotide == desired nucleotides

            def __choosehomo_or_wt__(targ_refalt,ansl,full_answer):#targrefalt, which ansl is it?
                                                                   #if ref matches targ[0] and targ[1] its WT+
                def __getansnuc__(an):                             #if ref matches targ[0] and alt matches targ[1]
                    annuc = an.split('-')                          #if ref matches none, its homovar
                    return(annuc[1:])

                if __getansnuc__(ansl)[0] != __split__(nucl)[0] and  __getansnuc__(ansl)[1] == __split__(nucl)[1]:#ref ne inc1 and both nucs are same
                    return full_answer['homohgvs']

                elif __getansnuc__(ansl)[0] == __split__(nucl)[0] and  __getansnuc__(ansl)[1] == __split__(nucl)[1]:#this is the ref alt case, ref = nuc1 and alt = nuc2
                    return full_answer['hethgvs']

                elif __getansnuc__(ansl)[0] == __split__(nucl)[1] and  __getansnuc__(ansl)[1] == __split__(nucl)[0]:#this is the alt ref case,
                    return full_answer['hethgvs']

                elif __getansnuc__(ansl)[0] == __split__(nucl)[0] and  __getansnuc__(ansl)[0] == __split__(nucl)[1]:
                    return full_answer['wthgvs']


            targ_refalt = __split__(nucl)#list containing the target ref  --   alt
            correct_hgvs = __choosehomo_or_wt__(targ_refalt,ansl,full_answer)#get the correct hgvs in the ansl

            def __fixhgvs__(hgvs,dval,nuc,desval):#straight split and reverse??yes.
                try:
                    hgvslist = hgvs.split(':')
                    newhgvs = ':'.join(hgvslist[1:]) + "||" + hgvslist[0]
                    return newhgvs

                except AttributeError:
                    print(str(dval) + "\t" + str(hgvs) + "\t" + str(nuc) + "\tUNABLE TO FIND RSID MAP" + "\t" + str(desval))

            fixedhgvs = __fixhgvs__(correct_hgvs,targ,nucl,desval)
            return(fixedhgvs)

        #---main control 4or4 i__match__rs___   ----

        for ansline in answer:#construct final returns two values as there are two matches, unfound because there is a mismatch,
                              #NEXT:: make specific for the alt value needed for multi alts, and reverse alts
            def __choose_correctans__(rsinc,rsans):#return the correct ALT
                                                   #find if the alt matches alt.
                altval = "NoneText"                #ref alt ONLY for the correct answer line
                def __get_REF_rsans__(ans):
                    refgetter = ans.split('-')

                    return refgetter[1]

                inc_nuc = __split__(rsinc)
                ref =__get_REF_rsans__(rsans)

                if ref in rsinc and inc_nuc[0] != inc_nuc[1]:#this is to test if this is het
                    #print("IN THE HETTER  " + str(ref) + "\t" + str(inc_nuc))

                    if ref == inc_nuc[0]:
                        altval = inc_nuc[1]

                    else:
                        altval = inc_nuc[0]

                elif ref in rsinc and inc_nuc[0] != inc_nuc[1] and ref == inc_nuc[1]:#ref is inc_nuc[1]
                    altval = inc_nuc[1]

                elif ref in rsinc and inc_nuc[0] != inc_nuc[1] and ref == inc_nuc[1]:#ref is inc_nuc[0]
                    altval = inc_nuc[0]

                elif ref not in rsinc:
                    altval = inc_nuc[0]

                else:
                    altval = ref

                return altval

            if rs_targ['rsID'] in ansline['rsid'] and __choose_correctans__(rs_targ['Diplotypes'],ansline['rsid']) in ansline['rsid']:#here, filter for correct

                correctcdot = __matchalt__(rs_targ,ansline['rsid'],desiredval['Diplotypes'],ansline,desiredval)
                #print(str(correctcdot) + " should be ok for correct dot")
                return correctcdot

    def __get_trans__(sym,trfi):#get trans, with read to new dict as sym for key
        trfi.seek(0)

        def __createdict__(trf):
            trnsdict = {}
            for line in trf:
                trnsdict[line['sym']] = line['trans']
            return trnsdict

        transcriptid = __createdict__(trans)[sym]
        return transcriptid

    if desiredval['rsID'] is not "":#match RSID to cdot, then to correct FULL HGVS.
        #variant_uuid = __match_rs__(desiredval['rsID'])#change, there are multiple
        variant_uuid = __match_rs__(desiredval)#change, there are multiple

    elif desiredval['Diplotypes'] is not "":#give justhe NM_ID for haplotypes here. also, set the NAT2. becareful of *N/*M appearing as *M/*N
        variant_uuid =  desiredval['Diplotypes'] +"||"+ __get_trans__(desiredval['Gene'],transf)

    else:#assuming its just the naked genetype
        variant_uuid =  "ALL||" + __get_trans__(desiredval['Gene'],transf)

    return variant_uuid


def printvar(met,vape,outfile):



    outfile.writerow([vape["Write Up ID"],vape['UUID'],vape["Drug Name"],met])

    #print(vape['UUID'])

#####----------------MAIN--------------####      #####----------------MAIN--------------####

mapdat = {}

#construct map for replacement table
for mapln in maps:

    valid_uuid = getVal(mapln,answer)
    mapdat[valid_uuid] = mapln['Metabolizer Status']

#map constructed, proceed.


newvapor.writerow(["Write Up ID","UUID","Drug Name","Metabolizer Status"])

for vape in vapors:# OK, well, get the uuid, then, decide if it is NAT2, or not NAT2.
                   #then decide if it is SNPs or not SNPs finally, get the HAPs.
                   #access the all if needed, but just make sure the Metabolizer Status is filled in.
                   #
                   #only write out the things that are in:
                   #    # error
                        #SNPs
                        #HAP2 snps

    def __getUUID__(vp):
        if "c.1865" in vp['UUID']:
            tmpid = vp['UUID'].split('+')
            ruuid = '+'.join(tmpid[0:2])
            return(ruuid)

        else:
            #print(vp['UUID'])
            ruuid = vp['UUID'].split('+')[0]
            #print(ruuid)
            return(ruuid)

    vuuid = __getUUID__(vape)#

    if 'NM_000015' in vuuid and "*" in vuuid:#get the digits try forward and backwards.
        g = re.search('(\*\d+)/(\*.+?)\|\|(.+)',vuuid)

        alsouuid = g.group(2) + "/" + g.group(1) + "||" + g.group(3)

        try:
            met_status = mapdat[alsouuid]#mapdat is the met status from the file.
            printvar(met_status,vape,newvapor)

        except KeyError:
            met_status = mapdat[vuuid]#mapdat is the met status from the file.
            printvar(met_status,vape,newvapor)


    elif ":LC" in vuuid:
        #newvuuid = re.sub(r':LC',"",vuuid)
        #met_status = mapdat[newvuuid]#mapdat is the met status from the file.
        printvar("This genotype could not be determined",vape,newvapor)

    elif ":Novel" in vuuid:
        #newvuuid = re.sub(r':Novel',"=",vuuid)
        #met_status = mapdat[newvuuid]#mapdat is the met status from the file.
        printvar("This genotype could not be determined",vape,newvapor)

    elif 'c.' in vuuid:

        try:
            met_status = mapdat[vuuid]#mapdat is the met status from the file.
            printvar(met_status,vape,newvapor)

        except KeyError:
            print("Unable to find Metabolizer Status for\t\""+vuuid+"\"\t"+str(vape['Drug Name'])+"\t"+str(vape['Evidence'])+"\t"+str(vape['Recommendation'])+"\t"+str(vape['Implication'])+"\t"+str(vape['Severity'])+"\t"+str(vape['Archived Implication Field']))



    else:
        print("the fuck is this, just make sure met status in not empty ")




















