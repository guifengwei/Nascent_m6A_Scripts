#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 27 Jul 2020 15:39:53

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
        for i, intron in enumerate(Introns):
            if bed.overlap(bed, intron):
                L = bed.overlapLength(intron)
                f = bed_in_Intron(Bed(line[0:6]), [intron])
                overlap_start, overlap_end = max(bed.start, intron.start), min(bed.end, intron.end)
                try:
                    random_start = random.randint(0, intron.length() - 1 - L) + intron.start
                    random_end  = random_start + L
                except:
                    while i>0:
                        X = Introns[i-1]
                        if X.length() >= L+1:
                            random_start = random.randint(0, X.length() - 1 - L) + X.start
                            random_end  = random_start + L
                            break
                        else:
                            i = i-1
                if f+1>0:
                    print( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(bed, intron, f, overlap_start, overlap_end, random_start, random_end, intron.length(), Introns[random_n].length()) )

if __name__ == "__main__":
    main()


