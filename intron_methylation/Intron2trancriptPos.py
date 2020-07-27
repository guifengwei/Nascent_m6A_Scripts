#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 27 Jul 2020 15:40:21

import sys, os, argparse, string, random
from UCSC_Class import Bed, genebed

def summitTobed(bedfile):
    '''
    summit 2 bed_interval
    '''
    summit2bed={}
    for line in open(bedfile, "r"):
        line = line.strip().split("\t")
        bed=Bed(line[0:6])
        summit2bed[bed.id] = bed
    return summit2bed

def bed_in_Intron(bed, Introns):
    '''return the relative position into individual intron, without cumulation '''
    for intron in Introns:
        if bed.overlap(bed, intron):
            f = 1.0*(bed.end - intron.start)/intron.length()
            if bed.strand == intron.strand and bed.strand == "+":
                return f
            elif bed.strand == intron.strand and bed.strand == "-":
                return 1-f
            else:
                return -1
        return -1

def main():
    ''' main scripts '''
    summit2bed=summitTobed("../m6AIP_Cfg1_Cfg2.narrowPeak.refined.bed")
    for line in open(sys.argv[1], 'r'):
        line = line.strip().split("\t")
        ### from the bed.id to bed interval
        try:
            bed=summit2bed[line[3]]
        except:
            sys.exit("# Something is wrong with the summit2bed")
        gene=genebed(line[6:18])
        Introns = gene.Introns()
        random_n = random.randint(0, len(Introns)-1)
        f1=bed_in_Intron(Bed(line[0:6]), [Bed(line[6:12])])
        random_start = random.randint(gene.start, gene.end)
        f2=bed_in_Intron(Bed([bed.chr, random_start, random_start+1, "None", 0, bed.strand]), [Bed(line[6:12])])
        print( "{0}\t{1}".format(f1,f2))

if __name__ == "__main__":
    main()


