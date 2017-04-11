#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/vbin/PyVGRes')
import vgr
import csv
from subprocess import call
import pybedtools
import random
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map;
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
parser = argparse.ArgumentParser(description="readthrough the omicia CSV, and if needed, reannotate the vcfs for the correct MED.res.files. now ncludes the RESULTS.txt library!! \n\n THIS is to recover the PREV files WITH NO One-OFF fix ")
parser.add_argument("--csv",help="the omicia csv file(s).")
parser.add_argument("--vcf",help="fullvcf file")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file, bgzipped and tabixed as -p vcf")
parser.add_argument("--sample",help="sample name")
parser.add_argument("--skip",help="list of items to skip")
parser.add_argument("--bed",help="the genome bed file to report on")

args = parser.parse_args()
vcf1fi = args.csv
vcffi = args.vcf
resfi = args.res
sample = args.sample
skipfi = args.skip
bedfi = args.bed

#special load for vcffi? (vcffi,"r")
omicia = csv.DictReader(open(vcf1fi,'r'))
vcf_full = vcf.Reader(open(vcffi,'r'))
#vcf_writer = vcf.Writer(open('redo.Merged.vcf', 'w'), vcf_full)
newres = vgr.Writer(open('NEW.' + sample + '.RESULTS.txt',"w"))
res = vgr.Reader(open(resfi,'r'))
reporter = open(sample + ".REPORT.txt","w")
skiplist = csv.DictReader(open(skipfi,"r"))


results = {}#stores the results lines;
omiciain = 0
recovered = 0
oneoffed = 0
skipct = 0
#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

def LoserRecover(ovcf,rsid):
    failreturn = vgr.Reader("10")#instanciate empty record
    newresults = vgr.Reader(open(rsid + '.COMPLETE.txt',"r"))
    success = 0
    return_value = ""
    for reslines in newresults:
        formattedlines = AddOmicia(ovcf,reslines)
        if formattedlines:
	    formattedlines
            success += 1
            return_value = formattedlines

    if success == 0:
        #print ovcf.ALT
        failreturn.CHROM = ovcf[0]
        failreturn.POS = str(ovcf[1])
        failreturn.REF = 'FABRICated'
        failreturn.ALT = 'VAR'
        failreturn.INFO = {}
        failreturn.INFO['VEP_EFFECT'] = "FABRIC_HARD_TO_MAP(" + ovcf[3].rstrip() + ")|(" + ovcf[4].rstrip() + ")"
        failreturn.INFO['FBRefAlleleCount'] = 0
        failreturn.INFO['FBReferenceAlleleQ'] = 0
        failreturn.INFO['QUAL'] = ovcf[4].rstrip()
        failreturn.INFO['RSID'] = ovcf[3].rstrip()
        return_value =  failreturn
    return return_value



def AddOmicia(vcf,results):
    if 'VEP_EFFECT' in results.INFO.keys():
        results.INFO['VEP_EFFECT'] = results.INFO['VEP_EFFECT'].rstrip() + "(" + vcf[3].rstrip() + ")|(" + vcf[4].rstrip() + ")"
    elif 'EFF_EFFECT' in results.INFO.keys():
        results.INFO['EFF_EFFECT'] = results.INFO['EFF_EFFECT'] + "(" + vcf[3].rstrip() + ")|(" + vcf[4].rstrip() + ")"
    else:
        results.INFO['VEP_EFFECT'] = "(" + vcf[3].rstrip() + ")|(" + vcf[4].rstrip() + ")"
    results.INFO['RSID'] = vcf[3].rstrip()
    results.INFO['QUAL'] = vcf[4].rstrip()
    return results

def LoserWrite(record,rsid,name):#prot:
    filename = str(rsid) + ".Merged.vcf"
    vcf_writer = vcf.Writer(open(filename, 'w'), vcf_full)
    vcf_writer.write_record(record)
    #print rsid

def LoserReRun(record,rsid,name):
  #  command = "/vbin/GoPipeRUN/goVarAnnotateAndID.sh" + rsid
    #call(["/vbin/GoPipeRUN/goVarAnnotateAndID.sh",rsid])
    effout = open(rsid + '.efd.vcf',"w")
    call(['java -jar /vbin/snpEff/snpEff.jar eff -i vcf -csvStats -hgvs hg19 ' + rsid + '.Merged.vcf'],shell=True,stdout=effout)

    call(['/vbin/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf -i ' + rsid + '.Merged.vcf -o ' +rsid + '.VEP.vcf --everything --species homo_sapiens --cache --refseq --offline --assembly GRCh37 --fasta /vbin/ref/hg19.nochr.fa'],shell=True)
    complete = open(rsid + '.COMPLETE.txt',"w")
    call(['/vbin/Perl/matchPathsAndMergeCallers.2.pl  --evs /ref/EVS_AF.vcf -p /ref/FULL.GENOME.lookUP.txt -eff ' + rsid + '.efd.vcf -vep ' + rsid + '.VEP.vcf -trn /ref/WG_good_genes.lst  -exac /ref/ExAC.r0.3.1.sites.af.vcf  -hgmd /ref/Homo_sapiens.HGMD.hg19.chr.vcf'],shell=True,stdout=complete)

def rsornone(rsid):
    if rsid:
        return rsid
    else:
        return "rnd" + str(random.randrange(10000,999999))

def tempFileWriter(tmp_handle,name):
    handle = open(name,'w')
    handle.write('\n'.join(tmp_handle))




#####----------------MAIN--------------####      #####----------------MAIN--------------####




product_bed = pybedtools.bedtool.BedTool(open(bedfi,"r"))
#convert the skiplist to bed here, allow for readthrough.

skipper = []# pybedtools.bedtool.BedTool()
for skiprow in skiplist:
    start = int(skiprow['pos']) - 1
    stop = int(skiprow['pos']) + 1
    skipper.append(skiprow['chr'] + "\t" + str(start)  + "\t" + str(stop))

tempFileWriter(skipper,'skip.tmp')

skip_vars_bed = pybedtools.bedtool.BedTool(open('skip.tmp','r'))
#convert the omicia csv to annotated bed here:
omiciastring = []
ocount = {}#For the purpose of making the omicia unique and removing chrM
for oline in omicia:
    keycheck = oline['chromosome'] + oline['end_on_chrom']
    if not keycheck in ocount.keys():
        ocount[oline['chromosome'] + oline['end_on_chrom']] = 1
        om_st = int(oline['end_on_chrom']) - 1
        om_sp = int(oline['end_on_chrom']) +1
        line =  oline['chromosome'] + '\t' + str(om_st) + '\t'+ str(om_sp) + "\t" + rsornone(oline['rs_id']) +  "\t" + rsornone(oline['quality'])
        omiciastring.append(line)

tempFileWriter(omiciastring,'omicia.tmp')
omiciabed = pybedtools.bedtool.BedTool(open('omicia.tmp','r'))
#omicia_in_product = omiciabed.intersect(product_bed)
omicia_in_product = omiciabed.intersect(product_bed)
inproductct = len(omicia_in_product)
omiciain = inproductct
omicia_in_product = omicia_in_product.subtract(skip_vars_bed)
skipct =  inproductct - int(len(omicia_in_product))#vars that were removed from load


for omicia in omicia_in_product:
    omicia_line = str(omicia).split('\t')
    line_for_med = res.fetch(omicia_line[0],int(omicia_line[1]),int(omicia_line[2]))#pulline fr RESULTS.txt
    resrecord = ""
    for record in line_for_med:
        resrecord = record
    if resrecord:#NEXT, add info to res line
        res_to_write = AddOmicia(omicia_line,resrecord)
        newres.write_record(res_to_write)
        recovered += 1
    else: #Send to loser bracket, this is where magic happens and rerun the vcf line.
        vcfelements = vcf_full.fetch(omicia_line[0],int(omicia_line[1]),int(omicia_line[2]))
        vcfregion = None
        for vcftypes in vcfelements:#pull out each one in area:
            vcfregion = vcftypes
        if vcfregion is not None:#matches exact, get into new file and redo
            losermatch = 1
            rsid_holder = rsornone(omicia_line[3])
            #print str(vcfregion)
            LoserWrite(vcfregion,rsid_holder,sample)
            LoserReRun(vcfregion,rsid_holder,sample)
            line_2_add = LoserRecover(omicia_line,rsid_holder)

            newres.write_record(line_2_add)
            recovered += 1
        else:
            print "CANT FIND VAR\t" + str(omicia)
    #write the executor for the missing record.
        #if losermatch == 0: #only enter if NO EXACT MATCH. start with full region


report = "omicia\trecovered\toneoff\tfilteredout\n"
reporter.write(report)
reported = str(omiciain) + "\t" + str(recovered) + "\t" + str(oneoffed) + "\t" + str(skipct) + "\n"
reporter.write(reported)
