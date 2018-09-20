#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy
from collections import defaultdict
#exporter:https://vapor.veritasgenetics.com/?q=admin/structure/views/view/yt_default_search_solr/edit
#but really:https://vapor.veritasgenetics.com/?q=yt-get-all-variant-links-export&eid=5&return-url=yt-get-all-variant-links-export
parser = argparse.ArgumentParser(description='validation of the PGX load. takes the ')
parser.add_argument("--pgx", help="The preprepped upload file",required=True)
parser.add_argument("--ans", help="The product answer bed file",required=True)
parser.add_argument("--trns", help="The gene 2 var id mapper ",required=True)


args = parser.parse_args()

mutfi = args.ans
muts = csv.DictReader(open(mutfi,'r'),delimiter='\t')
toloadfile = args.pgx
loadfile = open(toloadfile,'r')

genefi = args.trns#GENE -> trans
geneload = open(genefi,'read')


##-----------defs--------##

def varPrint(row,value):#
    print row + "\t" + value




def isItWT(a1,a2,t1,t2):#construct expected, and match to that
    isit = False        #t1 is always

    inc = [a1,a2]
    targ = [t1,t1]

    if sorted(inc) == sorted(targ):
        isit = True

    return isit



def isItHOM(a1,a2,t1,t2):#Determines if this rsid pair is right for this call
    isit = False         #t1 is ref, only one that is set really.

    inc = [a1,a2]
    targ = [t2,t2]

    if sorted(inc) == sorted(targ):
        isit = True

    return isit



def isItHet(a1,a2,t1,t2):#if it is het, makes it a shit ton easier
    hets = False

    inc = [a1,a2]
    targ = [t1,t2]

    #print str(inc) +  "\t" + str(targ)

    if sorted(inc) == sorted(targ):
        hets = True

    return hets



def varMapper(rss,mainln,ans,gns):#rss is the incoming parsed rsid
    correct_hgvs = None


    grab = re.search('([ATGC]{2})\|\|(\w+)',mainln[3])

    nmtrans = grab.group(2)
    alleles = list(grab.group(1))#contains incoming alleles.


    if re.search('NM_\d+',nmtrans ) is None:
        nmtrans = gns[nmtrans]

    #next trick, assign the TWO allele vals to homo, or WT or HET

    for rsval in ans[rss]:#need to check both tho.

        #print rsval +  " top of each rsvalue "
        (rtoss,targ1,targ2) =  rsval.split('-')

        #strategy: the isit loops check if the call is correct with TorF
        #and if the
        #
        #

        ishet = isItHet(alleles[0],alleles[1],targ1,targ2)
        iswt = isItWT(alleles[0],alleles[1],targ1,targ2)
        ishom = isItHOM(alleles[0],alleles[1],targ1,targ2)

 #       print str(rsval) + "\tHET:" +  str(ishet) + "\tHOM:" + str(ishom) + "\tWT:" +   str(iswt)


        if ishet is True:
            correct_hgvs = ans[rss][rsval]['hethgvs']

        elif ishom is True:
            correct_hgvs = ans[rss][rsval]['homohgvs']

        elif iswt is True:
            correct_hgvs = ans[rss][rsval]['wthgvs']


    return correct_hgvs








def tree(): return defaultdict(tree)

#------------main--------##ok, here, read he vapordump into a dict. then, for each uuid, check to see if all paids are present.


genes = {}

for g in geneload:
    gn = g.split('\t')
    genes[gn[1].strip()] = gn[0]


answer = tree()


for line in muts:#load by rsid and load by genes

    rsid = line['rsid'].split('-')[0]
    #print rsid + "\t" + line['rsid']

    answer[rsid][line['rsid']] = line
    trns = line['homohgvs'].split(':')[0]

    answer[trns][line['rsid']] = line



for loadr in loadfile:  #goals here:validate that line was added.
    load = loadr.strip()#print correct loading line (uuid that is in and valid)
                        #HOW TO access things NOT loaded??ok. plan for each line, if in hte lookup list,
                        # print last call as IN. if not, print OUT
    ld = load.split(',')

    rs = re.compile('(rs\d+)')

    if rs.search(ld[3]) is not None:
        rsid = rs.search(ld[3]).group(0)
        varmapped  = varMapper(rsid,ld,answer,genes)

        if varmapped is not None:
            print str(varmapped) +"\t" + str(ld) + "\t" +  rsid



