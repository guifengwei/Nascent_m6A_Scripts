#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 01 Jul 2020 10:39:01

import sys, os, argparse, string
from collections import Counter, defaultdict
import re
from UCSC_Class import Bed, genebed
from pybedtools import BedTool

def sortbylength(bed):
    return bed.length()

def sortbycoordinate(bed):
    return bed.start

def Cluster2A5SSorA3SS(Cluster):
    ''' 'clu_5077': {'chrX\t166458208\t166462315': 0.0195360774825754, 'chrX\t166458343\t166462315': -0.0195360774825754}
    '''
    for k1,v1 in Cluster.items():
        if len(v1) == 2:
            bed2list=[]
            for k2,v2 in v1.items():
                bed2list.append(Bed(k2))
            bed2list.sort(key=sortbylength)
            short, longer = bed2list
            alternative=""
            ss_type=""
            if short.start == longer.start:
                alternative= BedTool(short.chr+"\t"+str(short.end)+"\t"+str(longer.end)+"\t"+k1+"alternative"+"\t0\t"+short.strand, from_string=True)
            elif short.end == longer.end:
                alternative= BedTool(short.chr+"\t"+str(longer.start)+"\t"+str(short.start)+"\t"+k1+"alternative"+"\t0\t"+short.strand, from_string=True)
            else:
                pass
            X = alternative.intersect(m6A_bed)
            longer=longer.chr+"\t"+str(longer.start)+"\t"+str(longer.end)
            if len(X) >= 1:
                yield "A5SS_A3SS", "m6A", k1, v1[longer]
            else:
                yield "A5SS_A3SS", "nom6A", k1, v1[longer]

def Cluster2A5SSorA3SS_v2(Cluster):
    ''' work for the multiple
    '''
    for k1, v1 in Cluster.items():
        if len(v1)>2:
            bedlist=[]
            for k2, v2 in v1.items():
                bedlist.append(Bed(k2))
            bedlist.sort(key=sortbylength)
            starts = [i.start for i in bedlist]
            ends = [i.end for i in bedlist]
            if len(set(starts)) ==1 or len(set(ends))==1:
                short, longer=bedlist[0], bedlist[-1]
                if short.start == longer.start:
                    alternative= BedTool(short.chr+"\t"+str(short.end)+"\t"+str(longer.end)+"\t"+k1+"alternative"+"\t0\t"+short.strand, from_string=True)
                elif short.end == longer.end:
                    alternative= BedTool(short.chr+"\t"+str(longer.start)+"\t"+str(short.start)+"\t"+k1+"alternative"+"\t0\t"+short.strand, from_string=True)
                else:
                    pass
                longer=longer.chr+"\t"+str(longer.start)+"\t"+str(longer.end)
                m6A_string = m6AORnot(alternative)
                yield "A5SS_A3SS_"+str(len(bedlist)), m6A_string, k1, v1[longer]

def Cluster2ExonSkipping(Cluster):
    '''
    '''
    for k1, v1 in Cluster.items():
        if len(v1) == 3:
            bed3list=[]
            for k2,v2 in v1.items():
                bed3list.append(Bed(k2))
            bed3list.sort(key=sortbycoordinate)
            [bed1, bed2, bed3] = bed3list
            longer=bed2
            longer=longer.chr+"\t"+str(longer.start)+"\t"+str(longer.end)
            if bed1.start == bed2.start and bed2.end == bed3.end and bed1.end < bed3.start:
                alternative = BedTool("\t".join([bed1.chr, str(bed1.end), str(bed3.start), k1+"alternative", "0", bed1.strand]), from_string=True)
                if len(alternative.intersect(m6A_bed))>=1:
                    yield "SE", "m6A", k1, v1[longer]
                else:
                    yield "SE", "nom6A", k1, v1[longer]

def Cluster2MXE(Cluster):
    ''' '''
    for k1, v1 in Cluster.items():
        if len(v1) == 4:
            bed4list=[ Bed(k2) for k2, v2 in v1.items() ]
            bed4list.sort(key=sortbycoordinate)
            [bed1, bed2, bed3, bed4] = bed4list
            if bed1.start == bed2.start and bed2.end > bed3.start and bed3.start > bed1.end and bed3.end == bed4.end and bed4.start > bed2.end:
                yield v1

def m6AORnot(bed):
    ''' report m6A or nom6A, input are BedTool type
    '''
    X = bed.intersect(m6A_bed)
    if len(X)>=1:
        return "m6A"
    else:
        return "nom6A"

def Cluster2pIntronRetention(Cluster):
    ''''clu_5077': {'chrX\t166458208\t166462315': 0.0195360774825754, 'chrX\t166458343\t166462315': -0.0195360774825754}
    '''
    for k1,v1 in Cluster.items():
        if len(v1) >2:
            for k2,v2 in v1.items():
                Junction = BedTool(k2, from_string=True)
                Junction_Bed = Bed(k2.split("\t"))
                X = Junction.intersect(Genes, wb=True)
                m6A_string=[]
                splicing_type=[]
                if len(X)>=1:
                    for gene in X:
                        g = genebed(gene[3:])
                        for e in g.Exons():
                            if Junction_Bed.overlap(Junction_Bed, e):
                                x_l = Junction_Bed.overlapLength(e)
                                if x_l >=10 and x_l < e.length()-2: ## if ==, which means the entire exon is inside of the intron
                                    exon = BedTool(str(e), from_string=True)
                                    alternative=Junction.intersect(exon)
                                    splicing_type.append("pfIntronRetention")
                                    m6A_string.append(m6AORnot(alternative))
                if len(set(m6A_string))==1 and len(splicing_type)>0:
                    fo_pfIR.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("|".join(set(splicing_type)), "|".join(set(m6A_string)), k1, k2, v2))

def ParseClutsterFile(filename):
    ''' '''
    Cluster=defaultdict(dict)
    for line in open(filename, 'r'):
        if line.startswith("intron"): continue
        line = re.split(":|\t", line.strip())
        coordinate="\t".join(line[0:3])
        new_dict = {coordinate:float(line[-1])}
        Cluster[line[3]].update(new_dict)
    return Cluster

def main():
    ''' main scripts '''
    global m6A_bed, Genes, fo_pfIR
    Cluster = ParseClutsterFile(sys.argv[1])
    m6A_bed = BedTool("/usr/people/bioc1387/Project/ChrM6A-seq/Peaks/m6AIP_Cfg1_Cfg2.narrowPeak.refined.bed")
    print("## Loading GENCODE vM22 annotation", file=sys.stderr)
    Genes = BedTool("/usr/people/bioc1387/Project/mm10/Annotation/GENCODE_vM22/gencode.vM22.rmdup.genebed")
    print("## Loading GENCODE vM22 annotation Done", file=sys.stderr)
    fo_A5SS = open("Splicing_A5SSA3SS", "w")
    fo_ES = open("Splicing_ES", "w")
    fo_pfIR = open("Splicing_pfIR", "w")
    for t, m, clu, psi in Cluster2A5SSorA3SS_v2(Cluster):
        fo_A5SS.write("{0}\t{1}\t{2}\t{3}\n".format(t, m, clu, psi))
    for t, m, clu, psi in Cluster2A5SSorA3SS(Cluster):
        fo_A5SS.write("{0}\t{1}\t{2}\t{3}\n".format(t, m, clu, psi))
    for t, m, clu, psi in Cluster2ExonSkipping(Cluster):
        fo_ES.write("{0}\t{1}\t{2}\t{3}\n".format(t, m, clu, psi))
    #for x in Cluster2MXE(Cluster):
    #    print(x)
    Cluster2pIntronRetention(Cluster)
    fo_A5SS.close()
    fo_ES.close()
    fo_pfIR.close()

if __name__ == "__main__":
    main()

