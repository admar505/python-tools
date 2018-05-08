#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
from subprocess import call
import csv
import random
sys.path.append("/home/nucleo/lib/PyVGRes")
import vgr

parser = argparse.ArgumentParser(description="select infor from VCF")
parser.add_argument("--vcf",help="vcf file, bgzipped and tabixed as -p vcf",required=True)
parser.add_argument("--list",help="lost to search through, in CAP format",required=True)
parser.add_argument("--trans",help="allowed transcripts, format example NM_014946.3",required=True)
parser.add_argument("--tag",help="allowed transcripts, format example NM_014946.3",required=True,action='append')
parser.add_argument("--id",help="sample id, so that the bams can be called correctly",required=True)

args = parser.parse_args()
vcffi = args.vcf
varslst = args.list
trans = args.trans
tagarray = args.tag

grandoutfile = open(args.id + ".CAPRESULTS.tsv",'w')

vcf_full = vcf.Reader(open(vcffi,'r'))
necessary_vars = csv.DictReader(open(varslst,'r'),delimiter="\t")

#-------------------------------------here by DEFSgONS!!----------------------------------*
def recallVar(needline):
    bedtargetname = needline['Gene_Name'] +"."+ needline['Chr'] + "."+ needline['Pos'] + ".bed"
    bedtarget = open(bedtargetname,'w')


    pos = re.search('(\d+)\-(\d+)',needline['Pos'])

    if pos is not None:
        bedtarget.write(needline['Chr'] +"\t"+ str(int(pos.group(1)) - 1) +"\t"+ str(int(pos.group(2)) + 1))

    else:
        bedtarget.write(needline['Chr'] + "\t" + str(int(needline['Pos']) -1) +"\t" + str(int(needline['Pos']) + 1))

    bedtarget.close()

    vcfoutname = needline['Gene_Name'] +"."+ needline['Chr'] + "."+ needline['Pos'] + ".vcf"
    vcffile = open(vcfoutname,'w')

    call(['freebayes --targets ' +bedtargetname + ' -f /ref/hg19.fa --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic --ploidy 2 -b ' + args.id +"."+ needline['Chr']  + '.bam | /vbin/vcflib/vcfallelicprimitives --keep-info --keep-geno '],shell=True,stdout=vcffile)
    sys.stdout.flush()
    vcffile.close()

    vcfonevar = vcf.Reader(open(vcfoutname,'r'))


    for capvar in vcfonevar:

        grandoutfile.write(needline['Gene_Name'] + "\t" + needline['HGNC_ID'] + "\t"+ needline['Trans'] +"\t"+ needline['Chr'] +"\t"+ str(needline['Pos']) +"\t"+  needline['Required'] +"\t"+ str(capvar.QUAL) +"\t"+ str(capvar.POS) +"\t"+ capvar.samples[0]['GT'] + "\n")

def checkTags(complete,taggies):
    return_array = []

    for tag in taggies:
        try:
            return_array.append(complete[tag])
        except (AttributeError, KeyError) as e:
            return_array.append(tag + "_NOT_FOUND")

    return return_array

def mergeNWrite(vgrfile,caps,qval):
    completefile = vgr.Reader(open(vgrfile,'r'))

    for line in completefile:
        printarray = checkTags(line.INFO,tagarray)

        grandoutfile.write(caps['Gene_Name'] + "\t" + caps['HGNC_ID'] + "\t"+ caps['Trans'] +"\t"+ caps['Chr'] +"\t"+ caps['Pos'] +"\t"+  caps['Required'] +"\t"+ str(qval) + "\t" + "\t".join(printarray) + "\n")
        #+"\t"+ line.POS +"\t"+ line.INFO['EFF_HGVS'] + "\t"+ line.INFO['FBReferenceAlleleQ'] +"\t"+  line.INFO['FBGenoType'] +"\t"+ line.INFO['VAPOR_URL']



def kickAnnotations(vcfln,capln,vcf_4_header):
    returnvals = []

    #command = "/vbin/GoPipeRUN/goVarAnnotateAndID.sh" + rsid
    #call(["/vbin/GoPipeRUN/goVarAnnotateAndID.sh",rsid])
    rsid = capln['Chr'] + "."  + capln['Pos']
    vcf_write = vcf.Writer(open(rsid + '.Merged.vcf','w'),vcf_4_header)
    vcf_write.write_record(vcfln)
    sys.stdout.flush()
    vcf_write.close()

    effout = open(rsid + '.efd.vcf',"w")
    errout = open(rsid + '.stderrr',"w")

    call(['java -jar /vbin/snpEff/snpEff.jar eff -i vcf -csvStats -hgvs hg19 ' + rsid + '.Merged.vcf'],shell=True,stdout=effout,stderr=errout)
    sys.stdout.flush()
    effout.close()

    call(['/vbin/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf -i ' + rsid + '.Merged.vcf -o ' + rsid + '.VEP.vcf --dir /vbin/ensembl-tools-release-86/  --everything --species homo_sapiens --cache --refseq --offline --assembly GRCh37 --fasta /vbin/ref/hg19.nochr.fa'],shell=True)
    sys.stdout.flush()

    complete = open(rsid + '.COMPLETE.txt',"w")

    call(['/vbin/Perl/matchPathsAndMergeCallers.2.pl  --evs /ref/EVS_AF.vcf -p /ref/FULL.GENOME.lookUP.txt -eff ' + rsid + '.efd.vcf -vep ' + rsid + '.VEP.vcf -trn '+ trans + '  -exac /ref/ExAC.r0.3.1.sites.af.vcf  -hgmd /ref/Homo_sapiens.HGMD.hg19.chr.vcf'],shell=True,stdout=complete,stderr=errout)

    mergeNWrite(rsid + '.COMPLETE.txt',capln,vcfln.QUAL)



def getVCFLine(varln,vfile):

    pos = re.search('(\d+)\-(\d+)',varln['Pos'])

    if pos is not None:
        found_vcf = vfile.fetch(varln['Chr'],int(pos.group(1)) - 1,int(pos.group(2)) + 1)

    else:
        found_vcf = vfile.fetch(varln['Chr'],int(varln['Pos']) -1,int(varln['Pos']) + 1)

    return_line = None

    for captured in found_vcf:
        return_line = captured

    return return_line




#####----------------MAIN--------------####      #####----------------MAIN--------------####



for needed in necessary_vars:
    extracted = getVCFLine(needed,vcf_full)
    if extracted is not None:
        kickAnnotations(extracted,needed,vcf_full)


    else:
        recallVar(needed)







