#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/ec2-user/lib/PyVGRes')
import vgr
import csv
import random
parser = argparse.ArgumentParser(description="select infor from VCF")
parser.add_argument("--full",help="full results file, I think. best if Beee-gzipped and tabixed.")
parser.add_argument("--res",help="the final files, with VGR format.",required=True,action='append')
parser.add_argument("--primary",help="the final current new result file. with VGR format.",required=True)
parser.add_argument("--clinvarnew",help="the clinvar with header, will be csv",required=True)
parser.add_argument("--clinvarold",help="the clinvar with header. This is the old file, only one allowed at this time, will be csv",required=True)




args = parser.parse_args()

###<---file - parts-->###


fullfi = args.full          #large full size results file.
reslst = args.res           #all results files
primefi = args.primary      #the primary result file
clinnew = args.clinvarnew   #the latest clinvar file
clinold = args.clinvarold   #the older clinvar files

vcf_full = vgr.Reader(open(vcffi,'r'))


newclin = csv.






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























