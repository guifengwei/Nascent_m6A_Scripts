#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 27 Jul 2020 15:40:08

import sys, os, argparse, string
from UCSC_Class import Bed, genebed
from collections import defaultdict

def bed_in_bedlist(bed, bedlist):
    ''' 
    '''
    for i, exon in enumerate(bedlist):
        if bed.overlap(bed, exon):
            return exon.id
    return ""

def main():
    ''' main scripts '''
    bed2anno=defaultdict(list)
    bedid2bed={}
    for line in open(sys.argv[1], "r"):
        line = line.strip().split("\t")
        bed=Bed(line[0:6])
        if bed.id not in bedid2bed:
            bedid2bed[bed.id]=line[0:6]
        gene=genebed(line[6:18])
        s = bed_in_bedlist(bed, gene.Exons())
        if s=="":
            s = bed_in_bedlist(bed, gene.Introns())
        if s=="":
            sys.exit("### intersection is wrong")
        bed2anno[bed.id].append(s)
    for k,v in bed2anno.items():
        x="\t".join(bedid2bed[k])
        ss = []
        ss+=["E" for e in v if "Exon" in e]
        ss+=["I" for e in v if "Intron" in e]
        if "E" in set(ss) and "I" in set(ss):
            alternative="Alternative_ExonIntron"
        elif "E" in set(ss) and "I" not in set(ss):
            alternative="Constitutive_Exon"
        elif "I" in set(ss) and "E" not in set(ss):
            alternative="Constitutive_Intron"
        else:
            pass
        print("{0}\t{1}\t{2}\t{3}".format(x, "|".join(v), "".join(ss), alternative))

if __name__ == "__main__":
    main()

