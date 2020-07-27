#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 18 Jun 2020 18:56:15

import sys, os, argparse, string
from collections import defaultdict

def readfile(filename, sep="\t"):
    ''' '''
    for line in open(filename, "r"):
        if not line.startswith("#"):
            line = line.strip().split(sep)
            yield line

def main():
    ''' main scripts '''
    prefix="m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_MaxORF_LongestNcRNA.intron."
    ##########
    peak2ExonIntron_anno = {}
    for x in readfile("m6AIP_Cfg1_Cfg2.summits.Intersects_with_GENCODE_vM24_AllTranscripts.ExonIntron_anno"):
        peak2ExonIntron_anno[x[0]]=x[1]
    print("## Loading the ExonIntron_anno, Done!", file=sys.stderr)
    ########
    peak_anno2Phastcon=defaultdict(list)
    for x in readfile(prefix+"Overlap.bed.Phastcon.bed"):
        peak_anno2Phastcon[x[3]].append(x[6])
    for x in readfile(prefix+"matchControl.bed.Phastcon.bed"):
        peak_anno2Phastcon[x[3]].append(x[6])
    print("## Loading the Phastcon for Overlap and matchControl, Done", file=sys.stderr)
    #######
    peak_anno2GC=defaultdict(list)
    for x in readfile(prefix+"Overlap.bed.GC"):
        peak_anno2GC[x[0]].append(x[1])
    for x in readfile(prefix+"matchControl.bed.GC"):
        peak_anno2GC[x[0]].append(x[1])
    print("## Loading the GC for Overlap and matchControl, Done", file=sys.stderr)
    #####
    peak2m6A_intensity={}
    for x in readfile("../Signal_plot/m6APeaks_Cfg_intensity.txt"):
        peak2m6A_intensity[x[3]]= [x[8], x[9]]
    print("## Loading the m6A_intensity for all the peaks, Done", file=sys.stderr)
    #####
    print("chr\tstart\tend\tpeak_id\tscore\tstrand\tintron_chr\tintron_start\tintron_end\tintron_id\tintron_score\tintron_strand\trPos\tOverlap_start\tOverlap_end\tControl_start\tControl_end\tIntron_length\tControl_length\tExonORIntron\tOverlap_Phastcon\tControl_Phastcon\tOverlap_GC\tControl_GC\tClass\tm6A_Intensity")
    for x in readfile(prefix+"Overlap_and_matchControl2"):
        peak=x[3]
        peak_anno=x[3]+"|"+x[9]
        ExonIntron_anno=peak2ExonIntron_anno[peak]
        Phastcon = peak_anno2Phastcon[peak_anno]
        GC = peak_anno2GC[peak_anno]
        Intensity = peak2m6A_intensity[peak]
        print( "{0}\t{1}\t{2}\t{3}\t{4}".format("\t".join(x),  ExonIntron_anno, "\t".join(Phastcon), "\t".join(GC), "\t".join(Intensity)) )

if __name__ == "__main__":
    main()

