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

def getVal(desiredval,answer):

    variant_uuid = None

    def __match_rs__(rs_targ):
        #find rsid in answer:
        ansfi.seek(0)

        def __matchalt__():





        for ansline in answer:

            if rs_targ in ansline['rsid']:





                print(desiredval)
                print(ansline)



    variant_uuid = __match_rs__(desiredval['rsID'])


    return variant_uuid




def callNotFound(new_input):
    print "NEW_CALL_UNKNOWN\t" + new_input.strip()


#####----------------MAIN--------------####      #####----------------MAIN--------------####


mapdat = {}


#construct map for replacement table
for mapln in maps:

    getVal(mapln,answer)



#olddat = {}#stores the old calls. chr:pos

#for oln in olds:
    #print str(getVal(oln,'chrN',transdat))


#print "position\tnewcall\toldcall"

#for nln in news:
#    col = nln.split('\t')
#    loc = col[0] + ":" + col[1]
#
#    if loc in olddat and col[0] != 'chrN':
#        effnew = getVal(nln,'EFF_HGVS',transdat)
#        print loc +"\t" + effnew + "\t" + olddat[loc]
#
#    elif col[0] == 'chrN':
#        effnew = getVal(nln,'chrN',transdat)
#
#        if effnew in olddat and effnew is not None:
#            print loc +"\t" + str(effnew) + "\t" + olddat[effnew]

#        else:
#            callNotFound(nln)

#    else:
#        callNotFound(nln)




