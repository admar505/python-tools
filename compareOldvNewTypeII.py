#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv

parser = argparse.ArgumentParser(description="for TypeII RESULTS files, compare new with old. thanks")
parser.add_argument("--old",help="previous call file, bgzipped and tabixed as -p vcf",required=True)
parser.add_argument("--new",help="fresh call file",required=True)
parser.add_argument("--trns",help="file of transcripts, gene names, and positions.")
args = parser.parse_args()
oldfi = args.old
newfi = args.new
bedfi = args.trns

olds = open(oldfi,'r')
news = open(newfi,'r')
trnss =  csv.DictReader(open(bedfi,'r'),delimiter='\t')



#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

def getVal(arr,TAG,trn):
    value = None
    cols = arr.split('\t')

    if cols[1] == 'NULL':
        cols[1] = '1000'


    if cols[0] == 'chr2' and cols[1] == '234668879':
        cols[1] = '234668919'

    def __getNM__(ar):

        returncall = None
        validnm = False
        geteff = re.search('EFF_HGVS=(\S+)',ar)
        getnm = re.search('(NM_\d+)\:\S+',geteff.group(1))

        if getnm is not None:
            returncall = getnm.group(1)
            validnm = True
        else:
            returncall = geteff.group(1)

        return [validnm,returncall]


    if TAG == 'location':

        if cols[0] == 'chrN':
            nmid = __getNM__(arr)
            if nmid[0] is True:
                value = ":".join(trn[nmid[1]])
            else:
                value = "chrN:1000"
        else:
            value = cols[0] + ":" + cols[1]

    elif TAG == 'EFF_HGVS':
        getnm = re.search('(EFF_HGVS=(NM_\d+\:\S+))',arr)

        if getnm is not None:
            value = getnm.group(2)

        else:
            value = cols[4]

    return value




#####----------------MAIN--------------####      #####----------------MAIN--------------####

transdat = {}

for line in trnss:
    transdat[line['trans']] = [line['chr'],line['start']]


olddat = {}#stores the old calls. chr:pos


for oln in olds:
    olddat[getVal(oln,'location',transdat)] = getVal(oln,'EFF_HGVS',transdat)
    #print str(olddat)

print "position\tnewcall\toldcall"

for nln in news:
    col = nln.split('\t')
    loc = col[0] + ":" + col[1]
    if loc in olddat:
        effnew = getVal(nln,'EFF_HGVS',transdat)
        print loc +"\t" + effnew + "\t" + olddat[loc]


    else:
        print "NEW Call NOT FOUND:" +  nln.strip()















