#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy
import hgvs
from hgvs import *

import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.config
from hgvs import variantmapper

#from hgvs.variantmapper import  AssemblyMapper


##This is version 2, will now do the lookup in HGVS lands, using local HGVS
#    # after install    #    #
#export HGVS_SEQREPO_DIR=/usr/local/share/seqrepo/2018-11-26
#export UTA_DB_URL=postgresql://anonymous@localhost:15032/uta/uta_20170117
##

parser = argparse.ArgumentParser(description='given a list of hgvss\'/rsids from snpeff and VEP, see if the mapping is correct')
#parser.add_argument("--mut", help="Batch output",required=True)
parser.add_argument("--vg", help="the initial file, chr\tpos\tHGVS",required=True)
parser.add_argument("--rsid", help="if the file (--vg) is rsids, use this flag.",default = False)
args = parser.parse_args()

#mutfi = args.mut
vgrfi = args.vg
vgr = csv.DictReader(open(vgrfi,'r'),delimiter='\t')
#muts = open(mutfi,'r')


database = hgvs.dataproviders.uta.connect("postgresql://anonymous@localhost:15032/uta/uta_20170117")



##-----------defs--------##

def checkMatch(inc_hgvs,inc_chr,inc_pos,matchd):#newtranslated_pos,dictionary

    found = False

    if inc_hgvs in matchd:
        chk_val = "chr" + inc_chr + ":" + inc_pos

        if chk_val == matchd[inc_hgvs]:
            found =True

    return found




#------------main--------##

#blah, turn the hgvs into a dict

lookup = {}


for hg in vgr:

    isitgood = None

    varparse = hgvs.parser.Parser()
    vartomap = varparse.parse_hgvs_variant(hg['hgvs'])

    varmapper = hgvs.assemblymapper.AssemblyMapper(database, assembly_name='GRCh37', alt_aln_method='splign')


    mappedvar = varmapper.c_to_g(vartomap)

    print(mappedvar)




    #if isitgood == True:
    #    print("MATCH_FOUND\t" + resvalues[0] + "\tchr" + m.group(1) + "\t" + storelocale[storek[0]])

    #else:
    #    print("NOT_FOUND or NOT_RETURNED\t" + resln )





