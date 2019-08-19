#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv

parser = argparse.ArgumentParser(description="for TypeII RESULTS files, compare new with old. thanks")
parser.add_argument("--vapor",help="Current VAPOR dump",required=True)
parser.add_argument("--map",help="the file for translating the Metabolizer status",required=True)
parser.add_argument("--answer",help="Answer file, maybe. ")
args = parser.parse_args()
oldfi = args.vapor
newfi = args.map
ansfi = args.answer

mapfi = open(newfi,'r')
maps = csv.DictReader(mapfi,delimiter="\t")

ansfi = open(ansfi,'r')
answer = csv.DictReader(ansfi,delimiter="\t")

vaporfi =  open(oldfi,'r')
vapors =  csv.DictReader(vaporfi)
csvfile = oldfi + ".new.csv"
newvapor = csv.writer(open(csvfile, 'w'),  delimiter=',')


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
#


#-------------------------------------here by DEFSgONS!!----------------------------------*

def getVal(desiredval,answer):#the line being needed to map == desiredval, the full answerfile == answer

    variant_uuid = None

    def __match_rs__(rs_targ):#get correct rsid

        ansfi.seek(0)

        def __matchalt__(targ,ansl,nucl,full_answer,desval):#get the right hgvs nucleotide == desired nucleotides

            def __split__(nucl):#decide if snp or del
                if len(nucl) == 2:
                    return(list(nucl))

                else:
                    alts = nucl.split("/")#here, I am getting Del worked out

                    if 'D' in alts[0]:
                        return(alts[0].lower(),alts[1].lower())

                    elif 'D' in alts[1]:
                        return(alts[0],alts[1].lower())

                    else:
                        return(alts[0],alts[1])

            def __choosehomo_or_wt__(targ_refalt,ansl,full_answer):#targrefalt, which ansl is it?
                                                                   #if ref matches targ[0] and targ[1] its WT+
                def __getansnuc__(an):                             #if ref matches targ[0] and alt matches targ[1]
                    annuc = an.split('-')                          #if ref matches none, its homovar
                    return(annuc[1:])

                if __getansnuc__(ansl)[0] != nucl[0] and  __getansnuc__(ansl)[1] == nucl[1]:
                    return full_answer['homohgvs']

                elif __getansnuc__(ansl)[0] == nucl[0] and  __getansnuc__(ansl)[1] == nucl[1]:
                    return full_answer['hethgvs']


                elif __getansnuc__(ansl)[0] == nucl[0] and  __getansnuc__(ansl)[0] == nucl[1]:
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
            if rs_targ['rsID'] in ansline['rsid']:
                correctcdot = __matchalt__(rs_targ,ansline['rsid'],desiredval['Diplotypes'],ansline,desiredval)
                #print(str(correctcdot) + " should be ok for correct dot")
                return correctcdot

    if desiredval['rsID'] is not "":#match RSID to cdot, then to correct FULL HGVS.
        #variant_uuid = __match_rs__(desiredval['rsID'])#change, there are multiple
        variant_uuid = __match_rs__(desiredval)#change, there are multiple


    else:#give justhe NM_ID for haplotypes here. also, set the NAT2. becareful of *N/*M appearing as *M/*N
        #print(str(desiredval) + ' pass through')
        variant_uuid = 'pass through'


    return variant_uuid





#####----------------MAIN--------------####      #####----------------MAIN--------------####


mapdat = {}


#construct map for replacement table
for mapln in maps:

    print("INCOMING MAP LINE \t" + str(mapln))
    valid_uuid = getVal(mapln,answer)
    print(str(valid_uuid) + "\tASSIGNED UUID")



#olddat = {}#stores the old calls. chr:pos

#for oln in olds:
    #print str(getVal(oln,'chrN',transdat))


#print "position\tnewcall\toldcall"





