#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/PyVGRes')
import vgr
import csv
import random
parser = argparse.ArgumentParser(description="select infor from VCF")
parser.add_argument("--vgr",help="vcf file, bgzipped and tabixed as -p vcf",required=True)
parser.add_argument("--tag",help="tag to pull",required=True,action='append')

args = parser.parse_args()
vcffi = args.vgr
taglst = args.tag

vcf_full = vgr.Reader(open(vcffi,'r'))

#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def anyNone(rets):

    size = len(rets)
    none_ct = size  #set to size, reduce as nones go


    for tagkey in rets:
        checkfornone = re.match('.*?\[None\].*',str(rets[tagkey]))#LOGIC: if there is a none check for none will be NOT NONE??

        try:

            if checkfornone is  None:
                none_ct = none_ct - 1 #decided to use as --, as it is more sensible.

        except (AttributeError, IndexError) as e:
            dn = open(os.devnull,'w')

    final_return = None

    if none_ct < size:##indicates that no NONE vals were found in entirety.
        final_return = rets

    return final_return


def getTags(tags, varset):
    retval = {}

    for tagval in tags:

        if tagval in varset:
            #retval.append(tagval +  '=' + str(varset[tagval]))
            retval[tagval] = varset[tagval]

    #make a loop or def() that checks of at least one is not none.

    return_final = anyNone(retval)
    #print return_final
    return return_final


def makePLine(dct):
    returnvals = []
    for tags in dct.keys():
        returnvals.append(str(tags) + "=" + str(dct[tags]))

    return "\t".join(returnvals)


#####----------------MAIN--------------####      #####----------------MAIN--------------####

for line in vcf_full:
    good_tags = getTags(taglst,line.INFO)

    if good_tags is not None:#might need to manage each tag individually
        print line.CHROM + "\t" + str(line.POS) + "\t" + str(line.REF) + "\t" + str(line.ALT) + "\t" +   makePLine(good_tags)























