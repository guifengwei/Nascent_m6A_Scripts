#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 20 Jul 2020 12:29:02

import sys, os, argparse, string
from zbed import Bed
import math

def main():
    ''' main scripts '''
    Peaks = []
    for line in open(sys.argv[1], "r"):
        line = line.strip().split("\t")
        Peaks.append( Bed(line) )
    for i, bed in enumerate(Peaks):
        try:
            if bed.overlap(bed, Peaks[i-1]) or bed.overlap(bed, Peaks[i+1]):
                if bed.overlap(bed, Peaks[i+1]):
                    new_stop = math.floor( (int(bed.stop) + int(Peaks[i+1].start))/2 - 1 )
                else:
                    new_stop = bed.stop
                if bed.overlap(bed, Peaks[i-1]):
                    new_start = math.ceil( (int(bed.start) + int(Peaks[i-1].stop))/2 )
                else:
                    new_start = bed.start
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(bed.chr, int(new_start), int(new_stop), bed.id, int(bed.score), bed.strand))
            else:
                print(bed)
        except IndexError:
            print(bed)

if __name__ == "__main__":
    if len(sys.argv)==2:
        main()
    else:
        print("Usage: python refinePeak_neighbour.py MACS2_peaks.narrowPeak.bed")

