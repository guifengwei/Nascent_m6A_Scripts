#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 20 Jul 2020 12:27:41

import sys, os, argparse, string
from zbed import Bed 
import re

### fragmented length
f=125

def main():
    ''' main scripts '''
    Peak = []
    for line in open(sys.argv[1], "r"):
        line = line.strip().split("\t")
        Peak.append( Bed(line[0:4]) )
    ########################
    Strand={}
    Summit={}
    for line in open("m6AIP_abcam_summits.refined.bed", "r"):
        line = line.strip().split("\t")
        Strand[line[3]] = line[5]
        Summit[line[3]] = int(line[2])
    ###########################
    for i, bed in enumerate(Peak):
        bed.score = 0
        try:
            bed.strand = Strand[bed.id]
        except KeyError:
            continue
        if not re.findall("[a-z]$", bed.id):
            if bed.stop - bed.start <= f * 2:
                print(bed)
            else:
                bed.start = Summit[bed.id] - f
                bed.stop = Summit[bed.id] + f
                print(bed)
        else:
            if bed.start <= Summit[bed.id] - f:
                bed.start = Summit[bed.id] - f
            if bed.stop >= Summit[bed.id] + f:
                bed.stop = Summit[bed.id] + f
            print(bed)

if __name__ == "__main__":
    if len(sys.argv)==3:
        main()
    else:
        print("Usage: python refinePeak_fragLength.py MACS2.peaks.strand.narrowPeak MACS2.summits.strand.bed")


