#!/usr/bin/env python
import argparse
import subprocess
import re
import csv

parser = argparse.ArgumentParser(description='given a list of hgvss\'/rsids from snpeff and VEP, see if the mapping is correct')
parser.add_argument("--vg", help="the initial file, chr\tpos\tHGVS",required=True)
args = parser.parse_args()

vgrfi = args.vg
vgr =open(vgrfi,'r')

##-----------defs--------##

def retType(check_me):
    alt_type = None

    if 'indel' in check_me:
        alt_type = 'ins'

    elif 'ins' in check_me:
        alt_type = 'indel'

    elif 'del' in check_me:
        alt_type = 'del'

    elif 'dup' in check_me:
        alt_type = 'dup'

    elif '>' in check_me:
        alt_type = "snp"

    return alt_type


def checkMatch(cdot,gdotchange):#

    c_type = retType(cdot)
    g_type = retType(gdotchange)

    match = False#if False, the alteration type.
    if c_type == g_type:
        match = True

    found = [str(c_type), str(g_type), str(match)]

    return found


#------------main--------##
sum_all = {}

for hg in vgr:

    cols = re.split(r'\s+',hg)

    if cols[0]  not in sum_all:
        sum_all[cols[0]] = int(1)

    else:
        sum_all[cols[0]] = int(sum_all[cols[0]]) + 1




    if 'INCORRECT_POSITION' in cols[0]:


        if re.match('NC_',cols[4]):
            m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)(.*)',cols[4])
            distance = int(cols[3]) - int(m.group(3))

            changetype = checkMatch(cols[1],m.group(4))

            print(str(distance) + "\t" + "\t".join(changetype) + "\t" + cols[1] + "\t" + cols[2] +"\t" + cols[3] + "\t" + cols[4])

        else:
            print("Not working:"+ cols[1])



for sums in sum_all:
    print(sums + "\t" + str(sum_all[sums]))




