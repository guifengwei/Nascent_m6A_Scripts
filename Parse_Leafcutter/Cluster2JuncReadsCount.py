#!/usr/bin/python
# Guifeng Wei, <guifeng.wei@bioch.ox.ac.uk>
#-- coding: utf-8 --
#Last-modified: 26 Dec 2020 23:45:27

import pysam, argparse, sys, os
from collections import Counter,defaultdict
from UCSC_Class import Bed,genebed

def calculateJuncReads(bamfile, IntronBed):
    """
    extact the reads spaning or covering the splicing site,
    """
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    for intron in IntronBed:
        Junction_Reads = []
        counts = defaultdict(list)
        reads = bamfile.fetch( reference=intron.chr, start=intron.start, end=intron.end )
        for read in reads:
            if 'N' in read.cigarstring:
                blocks = read.get_blocks()
                Starts, Ends = zip(*blocks)
                if intron.start in Ends and intron.end in Starts:
                    counts["Exon_Exon"].append(read.query_name)
        intron.score = len(counts["Exon_Exon"])
        print(intron)

def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --cluster leafcutter_ds_effect_sizes.txt --bam bamfile', epilog="dependency: python3, UCSC_Class ")
    ############ GenePred format annotation
    p.add_argument('-c','--cluster', dest='cluster', metavar='Cluster',type=str, help="The Cluster file generated from LeafCutter: chr10:7801969:7802061:clu_6529")
    ########### RNA-seq alignment bamfile
    p.add_argument('-b', '--bam', dest="bamfile", metavar="Bam", type=str, help="The RNA-seq alignment bamfile")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def parseCluster(clusterfile):
    '''return the Bed of Junction, ie. intron.start and intron.end'''
    IntronBed = []
    with open(clusterfile, 'r') as f:
        for line in f:
            if not line.startswith("#") and not line.startswith("intron"):
                line = line.strip().split("\t")[0]
                junc = Bed(line.split(":"))
                junc.end -= 1
                IntronBed.append(junc)
    return IntronBed

def main():
    '''
    '''
    args = parse_argument()
    IntronBed = parseCluster(args.cluster)
    ################ Stranding the alignment
    #BAM_p, BAM_n = StrandingAlignment(args.bamfile, s=args.strand)
    print ("#### Starting to calculate the junction reads count ... ", file=sys.stderr)
    calculateJuncReads( bamfile=args.bamfile, IntronBed=IntronBed)
    print ("#### Done!", file=sys.stderr)

if __name__ == "__main__":
    main()

