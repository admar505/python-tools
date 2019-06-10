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

def checkMatch(inc_chr,inc_pos,mapd_chr,mapd_pos):#newtranslated_pos,dictionary

    found = False

    inc_val = inc_chr + ":" + inc_pos
    chk_val = "chr" + mapd_chr + ":" + mapd_pos

    if chk_val == inc_val:
    	found =True

    return found


#------------main--------##

#blah, turn the hgvs into a dict

varmapper = hgvs.assemblymapper.AssemblyMapper(database, assembly_name='GRCh37', alt_aln_method='splign')

for hg in vgr:

    storelocale = {}
    isitgood = None

    varparse = hgvs.parser.Parser()

    try:
        vartomap = varparse.parse_hgvs_variant(hg['hgvs'])
        mappedvar = varmapper.c_to_g(vartomap)

    except (hgvs.exceptions.HGVSDataNotAvailableError, hgvs.exceptions.HGVSInvalidIntervalError,hgvs.exceptions.HGVSParseError, hgvs.exceptions.HGVSInvalidVariantError) as e:

        print("NO_MAPPING\t" + str(hg['chr']) + "\t" + str(hg['pos']) + "\t" +  str(hg['hgvs']) )
        mappedvar = None


    if re.match('NC_',str(mappedvar)):
        m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)(.*)',str(mappedvar))
        storelocale[int(m.group(2))] = m.group(3) #loads the version, so that it can be sorted on.
        #print(m.group(3) + "\t" + m.group(2) + "\t" + m.group(1))
        #print("match " + m.group(2) )
        #print(str(len(storelocale.keys())) + "\t" +  str(storelocale.keys()))

    if len(storelocale.keys()) == 1:
        storek =  storelocale.keys()#this statement, will autosort??

        isitgood = checkMatch(hg['chr'],hg['pos'],m.group(1),m.group(3))



        if isitgood == True:
            print("MATCH_FOUND\t" + hg['hgvs'] + "\t" + hg['chr'] + "\t" + hg['pos']  + "\tchr" + m.group(1) + "\t" + storelocale[storek[0]])
            mappedvar = None
        else:
            print("INCORRECT_POSITION\t" + hg['hgvs'] + "\t" + hg['chr'] + "\t" + hg['pos']  + "\t" + str(mappedvar) )
            mappedvar = None




