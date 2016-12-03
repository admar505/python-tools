
import sys,os,re,fileinput,argparse
import vcf

parser = argparse.ArgumentParser(description="identify potential *5s")
parser.add_argument("--vcf",help="fullvcf file, or, just chr22 would be better but both work")
parser.add_argument("--t",help="ratio, as a float, ie,.56, in which to consider it to be *5+")

args = parser.parse_args()
vcffi = args.vcf
threshold = args.t
#special load for vcffi? (vcffi,"r")

vcf_full = vcf.Reader(open(vcffi,'r'))
hom_total = 0;
het_total = 0;
#use fetch to cut region:chr22:42522502-42540575
for cut_var in vcf_full.fetch('chr22',42522502,42540575):
    hom_total += cut_var.num_hom_alt
    het_total += cut_var.num_het

if het_total/(het_total + hom_total) < threshold:
    print "chr22\tNULL\tNULL\tFBGenoType=*5\tEFF_HGVS=*5\tRSID=NULL\tVAPOR_URL=http://vapor.veritasgenetics.com/?q=node/905420"

