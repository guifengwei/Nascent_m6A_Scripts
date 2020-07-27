#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# -- coding: utf-8 --

import sys,string
#from itertools import izip

''' convert genepred format to standard genebed format '''

class Bed:
    def __init__(self,x):
        if isinstance(x, str):
            x=x.rstrip("\n\r").split("\t")
        try:
            x[-1].rstrip("\n\r")
        except:
            pass
        self.chr=x[0].strip()
        self.start=int(x[1])
        if self.start<0:
            self.start=0
        self.end=int(x[2])
        try:
            self.id=x[3]
        except:
            self.id="NONAME"
        try:
            self.score=float(x[4])
        except:
            self.score=1
        try:
            self.strand=x[5]
        except:
            self.strand="."
        try:
            self.description="\t".join(x[6:])
        except:
            self.description=None
    def __str__(self):
        string=self.chr+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+str(self.id)+"\t"+str(self.score)
#       if(self.strand != "."):
#           string+="\t"+self.strand
        string+="\t"+str(self.strand)
        return string
    def length(self):
        return self.end-self.start
    def overlap(A,B):
        if A is None or B is None : return 0
        if(A.chr != B.chr) : return 0
        if (A.end <= B.start) : return 0
        if (B.end <= A.start) : return 0
        return 1
    overlap=staticmethod(overlap)
    def overlapLength(self, B):
        ''' return the overlap length between A and B '''
        if not B: return -1000000000
        if self.chr==B.chr:
            return self.length()+B.length()-max(self.end,B.end)+min(self.start,B.start)
        return -1000000000
    def merge(A,B):
        if A is None and B is None: return None
        if A and B: return Bed([A.chr,min(A.start,B.start),max(A.end,B.end),A.id+","+B.id,0,A.strand])
        if A: return A
        if B: return B
    merge=staticmethod(merge)
    def distance(self,bed):
        if(self.chr != bed.chr): return 10000000000
        if(Bed.overlap(self,bed)): return 0
        a=[abs(self.start-bed.end),abs(self.end-bed.end),abs(self.start-bed.start),abs(self.end-bed.start)]
        i=a[0]
        for j in a[1:]:
            if i>j:
                i=j
        return i
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.end,other.end) or cmp(self.strand,other.strand)
    def upstream(self,bp):
        '''return the $bp bp upstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_up"+str(bp)
        if(self.strand=="+"):
            start=self.start-bp
            end=self.start
        else:
            start=self.end
            end=self.end+bp
        if (start<0):start=0
        x=[chr,start,end,id,0,strand]
        return Bed(x)
    def downstream(self,bp):
        '''return the $bp bp downstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_down"+str(bp)
        if(self.strand=="+"):
            start=self.end
            end=self.end+bp
        else:
            start=self.start-bp
            end=self.start
        x=[chr,start,end,id,0,strand]
        return Bed(x)
    def upstreamextend(self,bp):
        '''return the bed itself and $bp bp upstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_upextend"+str(bp)
        if(self.strand!="-"):
            start=self.start-bp
            end=self.end
        else:
            start=self.start
            end=self.end+bp
        x=[chr,start,end,id,0,strand]
        return Bed(x)
    def downstreamextend(self,bp):
        '''return the bed itself and $bp bp downstream Bed Class Object'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_downextend"+str(bp)
        if(self.strand!="-"):
            start=self.start
            end=self.end+bp
        else:
            start=self.start-bp
            end=self.end
        x=[chr,start,end,id,0,strand]
        return Bed(x)
    def updownextend(self,bp):
        '''return the bed itself and $bp bp flanking regions Bed Class Object
'''
        chr=self.chr
        strand=self.strand
        id=self.id+"_updownextend"+str(bp)
        start=self.start-bp
        end=self.end+bp
        x=[chr,start,end,id,0,strand]
        return Bed(x)
    def strand_cmp(self,bed):
        if bed.strand == ".":
            return "."
        if self.strand == bed.strand:
            return "+"
        else:
            return "-"
    def contain(self,other):
        ''' Return 1 if self_bed contains other_bed '''
        if self.start <= other.start < other.end <= self.end:
            return True
        else:
            return False
    def LOCATE(self,B):
        ''' return the location relationship between B and self: IN, OVERLAP, INCLUDING, OUT, ANTISENSE '''
        if self.chr == B.chr:    
            if self.strand == B.strand:
                # B in self
                if self.start < B.start < B.end <= self.end or \
                            self.start <= B.start < B.end < self.end:
                    return "IN"
                # B ovelap with self
                if self.start <= B.start <= self.end < B.end or \
                            B.start < self.start <= B.end <= self.end:
                    return "OVERLAP"
                # B out of self
                if self.start > B.end or self.end < B.start:
                    return "OUT"
                # B include self
                if B.start <= self.start <= self.end <= B.end:
                    return "INCLUDING"
            else:
                if self.start <= B.start <= self.end or \
                            self.start <= B.end <= B.end:
                    return "ANTISENSE"
                else:
                    return "OUT"
        else:
            return "OUT"

class genepred():
    ''' genepred'''
    '''0-6: name(id), chrom, strand, txStart, txEnd, cdsStart, cdsEnd, \
       7-11: exonCount, exonStarts, exonEnds, score, name2, other(desc) '''
    def __init__(self,x):
        try:
            self.bin=int(x[0])
            if(self.bin < 10000):
                x=x[1:]
        except:
            pass
        self.id=x[0]
        self.chr=x[1]
        self.strand=x[2]
        if(self.strand == 1):
            self.strand="+"
        elif(self.strand == -1):
            self.strand="-"
        elif(self.strand == 0):
            self.strand="+"
        self.start=int(x[3])
        self.end=int(x[4])
        self.cds_start=int(x[5])
        self.cds_end=int(x[6])
        self.exon_count=int(x[7])
        self.exon_starts=x[8].split(",")
        self.exon_ends=x[9].split(",")
        for i in range(self.exon_count):
            try:
                self.exon_starts[i]=int(self.exon_starts[i])
                self.exon_ends[i]=int(self.exon_ends[i])
            except:
                pass
        self.exonSizes   = [self.exon_ends[i] - self.exon_starts[i] for i in range(self.exon_count) ]
        self.exonStartsR = [self.exon_starts[i] -self.exon_starts[0] for i in range(self.exon_count)]
        self.score=0
        try:
            self.protein_id=x[11]
            self.name2=x[11]
        except:
            self.name2="None"
        try:
            self.align_id=x[12]
        except:
            pass
    def __str__(self):
        return "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s,\t%s,\t%s" % (self.id,self.chr,self.strand,self.start,self.end,self.cds_start,self.cds_end,self.exon_count,",".join([str(self.exon_starts[j]) for j in range(self.exon_count)]),",".join([str(self.exon_ends[j]) for j in range(self.exon_count)]),self.name2)

class genebed():
    '''genebed '''
    ''' 0-5: chrom, chromstart, chromend, name, score, strand,\
        6-11: cdsStart, cdsEnd, itemRGB, exonCount, exonSizes, exonStarts '''
    def __init__(self,x):
        try:
            self.bin=int(x[0])
            if(self.bin<10000):
                pass
                #x=x[1:]
        except:
            pass
        self.chr = x[0].rstrip()
        self.start = int(x[1])
        self.end = int(x[2])
        self.id  = x[3]
        self.score = x[4]
        self.strand = x[5]
        self.cdsStart = int(x[6])
        self.cdsEnd = int(x[7])
        self.itemRGB= x[8]
        self.exonCount = int(x[9])
        self.exonSizes = x[10].strip(",").split(",")
        self.exonSizes = list( map(int, self.exonSizes))
        # coordinates relative to the transcript start(start)  self.exonStarts
        self.exonStarts = x[11].strip(",").split(",")
        self.exonStarts = list( map(int, self.exonStarts) )
        # coordinates in the chromosome
        # self.exon_starts self.exon_stops
        self.exon_starts = [each + self.start for each in self.exonStarts]
        self.exon_stops=[]
        for k,v in zip(self.exonStarts, self.exonSizes):
            self.exon_stops.append(self.start + k + v)
    def __str__(self):
       return "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s,\t%s," % (self.chr, self.start, self.end, self.id, self.score, self.strand, self.cdsStart, self.cdsEnd, self.itemRGB, self.exonCount, ','.join(map(str,self.exonSizes)), ','.join(map(str, self.exonStarts)) )
    def _exon(self,i):
        '''internal fucntion to call the exon position'''
        if i > self.exonCount:
            return None
        if self.strand=="+" :
            return self.exon_starts[i-1],self.exon_stops[i-1]
        if self.strand=="-":
            n=self.exonCount
            return self.exon_starts[n-i],self.exon_stops[n-i]
    def getExon(self,i):
        '''get exon and return Bed Class of exon, the first exon is Exon_1'''
        if(i > self.exonCount or i < 1):return None
        start,stop=self._exon(i)
        id=self.id+"_"+"Exon_"+str(i)
        x=[self.chr,start,stop,id,0,self.strand]
        return Bed(x)
    def Exons(self):
        '''return a list of Bed classes of exon'''
        a=[]
        for i in range(self.exonCount):
            a.append(self.getExon(i+1))
        return a
    def mRNA_length(self):
        '''sum of the length of exons'''
        s=0
        for i in self.Exons():
            s=s+i.length()
        return s
    def _intron(self,i):
        if i > self.exonCount-1:
            return None
        if self.strand=="+":
            return self.exon_stops[i-1],self.exon_starts[i]
        if self.strand=="-":
            n=self.exonCount
            return self.exon_stops[n-i-1],self.exon_starts[n-i]
    def getIntron(self,i):
        '''get exon and return Bed Class of intron, the first intron is Intron_1'''
        if(i>self.exonCount-1 or i< 1): return None
        start,stop=self._intron(i)
        id=self.id+"_"+"Intron_"+str(i)
        x=[self.chr,start,stop,id,0,self.strand]
        return Bed(x)
    def Introns(self):
        '''return a list of Bed classes of exon'''
        a=[]
        for i in range(self.exonCount-1):
            a.append(self.getIntron(i+1))
        return a
    def UTR5(self):
        ''' return the genebed class as UTR5 of some genes are also located in exon2 '''
        if(self.strand == "+"):
            if(self.cdsStart == self.start):
                return None
            return self.slice(self.start, self.cdsStart)
        if(self.strand == "-"):
            if(self.cdsEnd == self.end):
                return None
            return self.slice(self.cdsEnd, self.end) 
    def UTR3(self):
        ''' return the genebed class as UTR5 '''
        if(self.strand == "-"):
            if(self.cdsStart == self.start):
                return None
            return self.slice(self.start, self.cdsStart)
        if(self.strand == "+"):
            if(self.cdsEnd == self.end):
                return None
            return self.slice(self.cdsEnd, self.end)
    def CDS(self):
        ''' return the genebed class '''
        if(self.cdsStart == self.cdsEnd):
            return None
        return self.slice(self.cdsStart, self.cdsEnd)
    def region_cmp(self,bed):
        '''report Bed Object overlap with which exon and intron or don't overlap with gene'''
        if not Bed.overlap(self,bed):
            r="Non-overlap"
            return r
        r="Overlap: "
        for i in range(self.exonCount):
            start,stop=self._exon(i+1)
            if(bed.start<stop and start<bed.stop):
                r+="Exon_"+str(i+1)+","
        for i in range(self.exonCount-1):
            start,stop=self._intron(i+1)
            if(bed.start < stop and start < bed.stop):
                r+="Intron_"+str(i+1)+","
        return r
    def promoter(self,bp=1000):
        ''' Return the promoter region of the gene(Bed format) '''
        if not bp:
            bp = self.bp
        self.adjust_structure()
        if self.strand == '+':
            return Bed([self.chr,self.start-bp, self.start+bp, 'promoter_'+str(bp)+'_of_'+str(self.id),0,'+'])
        else:
           return Bed([self.chr,self.stop-bp, self.stop+bp, 'promoter_'+str(bp)+'_of_'+str(self.id),0,'-'])
    def slice(self, a, b):
        ''' get the genebed slice from coordinate a to coordinate b'''
        if self.start > a or self.end < b: return None
        x_start = a
        x_end = b
        x_cdsStart = a
        x_cdsEnd = b
        x_strand = self.strand
        #   x_cdsStart and x_cdsEnd
        if a<= self.cdsStart <= b:
            x_cdsStart = self.cdsStart
        elif self.cdsStart > b:
            x_cdsStart = b
        else:
            pass
        if a<= self.cdsEnd <= b:
            x_cdsEnd = self.cdsEnd
        elif self.cdsEnd < a:
            x_cdsEnd = a
        else:
            pass
        ##########################
        x_Exons=[]
        for k,v in zip(self.exon_starts, self.exon_stops):
            if a>v:
                continue
            elif k<=a<=b<=v:
                x_Exons.append([a, b])
                break
            elif k<=a and b>v:
                x_Exons.append([a, v])
                continue
            elif a<k and k<=b<=v:
                x_Exons.append([k, b])
                break
            elif a<k and b>v:
                x_Exons.append([k, v])
                continue
            elif b<k:
                break
        x_exonCount = len(x_Exons)
        x_exonSizes=""
        x_exonStarts=""
        for e in x_Exons:
            x_exonSizes += str( int(e[1]) - int(e[0]) )+","
            x_exonStarts+= str( int(e[0]) - a) + ","
        return genebed([self.chr, x_start,x_end, self.id+":"+str(a)+"_"+str(b), 0, x_strand, x_cdsStart,x_cdsEnd, self.itemRGB, x_exonCount,x_exonSizes,x_exonStarts])


