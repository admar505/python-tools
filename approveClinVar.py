#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/ec2-user/lib/PyVGRes')
import vgr
import vcf
import csv
import random
parser = argparse.ArgumentParser(description="select infor from VCF")
parser.add_argument("--prev",help="the prev result file, VGR format.",required=True)
parser.add_argument("--primary",help="the final current new result file. with VGR format.",required=True)
parser.add_argument("--acmg",help="the final current new result file. with VGR format.",required=True)
parser.add_argument("--clinvarnew",help="the clinvar with header, will be csv",required=True)
parser.add_argument("--clinvarold",help="the clinvar with header. This is the old file, only one allowed at this time, will be csv",required=True)




args = parser.parse_args()

###<---file - parts-->###


prevfi = args.prev          #all results files
primefi = args.primary      #the primary result file
clinnew = args.clinvarnew   #the latest clinvar file
clinold = args.clinvarold   #the older clinvar files
acmgfi = args.acmg


newclin = open(clinnew,'r')
oldclin = open(clinold,'r')

allres = vgr.Reader(open(prevfi,'r')) 
prime = vgr.Reader(open(primefi,'r')) 
acmg = open(acmgfi,'r')


#parse results in a map or dict, or what??
csv.field_size_limit(sys.maxsize)

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


def ClinSeek(clnvr,index):
    
    clnvr.seek(0)#sadly, I have to always reset to zero.

    def positiongather(chrom,start,stop,files):
        foundat = []
        for line in files:
            if str(line['#Chr']) == str(chrom) and int(line['Start']) >= int(start) and int(line['End']) <= int(stop):
                foundat.append(line['Alt'] + ":" + line['CLNSIG'] + ":" + line['hgmdVC'] + ":" + line['Start'] + ":" + line['End'])

        return(foundat)


    cln = csv.DictReader(clnvr,delimiter='\t')  #get the file


    chrom =  str(index.CHROM).split('r')[1:]    #get the tacos
    getcln = positiongather(chrom[0],int(index.POS) - 2, int(index.POS) +2,cln)
    
    return(getcln)

def measureindel(samp):#pull the ref and alt and return int(size of both) added? not sure yet.
    alt = len(samp.ALT)
    ref = len(samp.REF)
   
    return(int(alt) + int(ref)) 

def getposses(rec_ln):##just make returnable lines
    retln = str(rec_ln.POS) + "\t" +  str(rec_ln.ALT)
    return(retln) 


def getOld(index,sample):##$newsample, $oldsample   purpose, find the closest if not exact is what I am thinking.
                         #also, second is within the section of the edit, I am mainly thinking offset for the shift in indels.
    getres = []##list of positions found, one per                     
    #find locale, if not null, return infos.
    try:
        matchedlines = sample.fetch(index.CHROM,int(index.POS) - 1, int(index.POS) + 1)
    
        for line in matchedlines:
            getres.append(line)

    except ValueError:
        next 

        

    #check if it worked at all at all.

    if str(getres) == "[]":##OK here we go. IF the above didnt catch, then check if it is in the interval of the indel
        try:
            newwidths = measureindel(index)
            rescues = sample.fetch(index.CHROM,(int(index.POS) - newwidths), (int(index.POS) + newwidths))
    
            for rescue in rescues:
                getres.append(rescue)

        except ValueError:
            
            next

    return(getres)


def dummyfill(vals):

    valret = "."
    if vals:
      valret = vals

    return(valret)

def acmgGet(indxr,acmgf):

    getallacmg = []
    acmgf.seek(0)#sadly, I have to always reset to zero.
    acmgprepped = csv.DictReader(acmgf,delimiter='\t',fieldnames=["#Chr","Start","End","ref","Alt","type",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".",".","."])  #get the file

    def __positiongather__(chrom,start,stop,files):
        foundat = []
        for line in files:
           
            if str(line['#Chr']) == str(chrom) and int(line['Start']) >= int(start) and int(line['End']) <= int(stop):
                foundat.append(line['Alt'] + ":" +  line['Start'] + ":" + line['End'])

        return(foundat)





    tmpacmg = __positiongather__(indxr.CHROM,int(indxr.POS) -1, int(indxr.POS) +1,acmgprepped)
    for each in tmpacmg:
        getallacmg.append(each)

    return(','.join(getallacmg))



#####----------------MAIN--------------####      #####----------------MAIN--------------####

indexfoundinold = {}##if we foulnd this index in the old file, then record it here, we will go through the list and declare what lines were not found.
print("chrom\tref\tpos\tnew_alt\told_pos\told_alt\tNEW_Clinvar:ALT:clinpath:HGMD:CLNSTART:CLNSTOP\tNEW_Clinvar:ALT:clinpath:HGMD:CLNSTART:CLNSTOP\tACMG\tQUAL")
for line in prime:
    
    oldinfo = getOld(line,allres)
    newpath = dummyfill(ClinSeek(newclin,line)) 
    oldpath = dummyfill(ClinSeek(oldclin,line))
    acmgvals = dummyfill(acmgGet(line,acmg))
    
    foundkey = line.CHROM + ":" + line.POS    


    if str(oldinfo) == "[]":
        prl = getposses(line) 
        print(str(line.CHROM) +"\t" + str(line.REF) +"\t"+ str(prl) + "\tNO OLD\tMATCH\t" + ','.join((newpath))  + "\t" + ','.join((oldpath)) + "\t" +str(acmgvals) + "\t" + str(line.INFO['QUAL']))

    else:

        for info in oldinfo:
            prln = getposses(info)
            mainln = getposses(line)

            indexfoundinold[str(foundkey)] = foundkey
            
            print(str(line.CHROM) +"\t" + str(line.REF) + "\t" + str(mainln) + "\t" + str(prln) +"\t"+  ','.join((newpath))  + "\t" + ','.join((oldpath)) + "\t" + str(acmgvals) +"\t" + str(line.INFO['QUAL']))


print("BELOW\t are missing in the new run,\tthat were in old") 

oldres = vgr.Reader(open(prevfi,'r')) 


for clnup in oldres:#send this through, pick up what didnt catch, I think it is almost exactly like prev so may hafta figure elegant way.
    
    checkey = clnup.CHROM + ":" + clnup.POS

    if str(checkey) not in indexfoundinold:

        cnewpath = dummyfill(ClinSeek(newclin,clnup)) 
        coldpath = dummyfill(ClinSeek(oldclin,clnup))
        cacmgvals = dummyfill(acmgGet(line,acmg))
        

        print(str(clnup.CHROM) +"\t" + str(clnup.REF) + "\t" +str(clnup.POS)+ "\t" + str(clnup.ALT) + "\tMISSED IN\tNEW\t"+  ','.join((cnewpath))  + "\t" + ','.join((coldpath)) + "\t" + str(cacmgvals) +"\t"+ str(clnup.INFO['QUAL']))
