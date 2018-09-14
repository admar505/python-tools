#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy
#exporter:https://vapor.veritasgenetics.com/?q=admin/structure/views/view/yt_default_search_solr/edit
#but really:https://vapor.veritasgenetics.com/?q=yt-get-all-variant-links-export&eid=5&return-url=yt-get-all-variant-links-export
parser = argparse.ArgumentParser(description='validation of the PGX load. takes the ')
parser.add_argument("--pgx", help="the list from the loading of hte vapor",required=True)
parser.add_argument("--vap", help="the paid --> uuid list out of vapor",required=True)
args = parser.parse_args()

mutfi = args.vap
muts = open(mutfi,'r')
toloadfile = args.pgx
loadfile = open(toloadfile,'r')

##-----------defs--------##

def varPrint(row,value):#
    print row + "\t" + value


vapor={}

def addval(key,paid,vap):

    if key in vap:
        vap[key].append(paid)

    else:
        vap[key] = [paid]

def check(uuid,vap):#check to see if uuid format  was loaded
    in_loaded = False#could be called from in checkPA, but, didnt think of it soon enough

    if uuid in vap:
        in_loaded = True

    return in_loaded


def checkPA(valid,uuid,pa_id,vapor):#checks first if this uuid is valid, then, if the pa is valid. if not, returns False

    if valid is True:
        valid = False

        if pa_id in vapor[uuid]:
            valid = True

    return valid



#------------main--------##ok, here, read he vapordump into a dict. then, for each uuid, check to see if all paids are present.


for liner in muts:#should I load stuff in backwards? or forwards? lets try
    line = liner.strip()
    ln = line.split(',')

    addval(ln[0].strip(),ln[1].strip(),vapor)

for loadr in loadfile:  #goals here:validate that line was added.
    load = loadr.strip()#print correct loading line (uuid that is in and valid)
                        #HOW TO access things NOT loaded??ok. plan for each line, if in hte lookup list,
                        # print last call as IN. if not, print OUT
    ld = load.split(',')

    chkugt = re.compile('(c.\S+)(\[\d+\])\/(\[\d+\])(\|\|.+)')
    chkstar = re.compile('(\*\w+)\s?\/(\*?\w+)(\|\|.+)')

    ugt = chkugt.search(ld[3])
    star = chkstar.search(ld[3])

    if ugt is not None:#to be sure, I am checking the reverse val,
                       #this will help in future translations.
        revval = ugt.group(1) + ugt.group(3) + "/" +  ugt.group(2) + ugt.group(4)
        #check if forward val is in
        forward = check(ld[3],vapor)
        reverse = check(revval,vapor)

        uuid = ld[3]
        if forward is False:
            uuid = revval

        paforward = checkPA(forward,ld[3],ld[2],vapor)#ok, trick is if valid, then print line as yes,
        paback = checkPA(reverse,revval,ld[2],vapor)#but need to be able to read if both are false, to
                                                    #print no, not valid. but, also, print the valid
                                                    #lookup key, so, if forward and rev are false, nothing is good
        if forward is False and reverse is False:   #if paforward and paback are false, and  forward or rev is true,
            varPrint(load,"UUID_NOT_FOUND")         #then pa is not in, but uuid is.

        elif paforward is True:
            varPrint(load,"PA_ID_LOADED")

        elif paback is True:
            revline = ld[0] + "," + ld[1] + "," + ld[2] + "," + uuid + "," + ld[4] +","+ ld[5]
            varPrint(revline,"PA_ID_LOADED")

        elif paback is False and paforward is False:
            revline = ld[0] + "," + ld[1] + "," + ld[2] + "," + uuid + "," + ld[4] +","+ ld[5]
            varPrint(revline,"PA_NOT_LOADED")


    elif star is not None:
        revval = star.group(2) + "/" +  star.group(1) + star.group(3)

        forward = check(ld[3],vapor)
        rev = check(revval,vapor)

        paforward = checkPA(forward,ld[3],ld[2],vapor)
        paback = checkPA(rev,revval,ld[2],vapor)

        uuid = ld[3]
        if forward is False:
            uuid = revval

        if forward is False and rev is False:
            varPrint(load,"UUID_NOT_FOUND")

        elif paforward is True:
            varPrint(load,"PA_ID_LOADED")

        elif paback is True:
            revline = ld[0] + "," + ld[1] + "," + ld[2] + "," + uuid + "," + ld[4] +","+ ld[5]
            varPrint(revline,"PA_ID_LOADED")

        elif paback is False and paforward is False:
            revline = ld[0] + "," + ld[1] + "," + ld[2] + "," + uuid + "," + ld[4] +","+ ld[5]
            varPrint(revline,"PA_NOT_LOADED")

    else:
      print load.strip()












