#!/python

# Copyright (C) 2015 University of California, Los Angeles (UCLA)
# Shayna R. Stein, Emad Bahrami-Samani, Yi Xing
#
# Authors: Shayna R. Stein, Emad Bahrami-Samani, Yi Xing
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.

# standard python imports
import sys,re,os,shutil,gzip
import re,logging,time,datetime,commands,argparse,random
from collections import defaultdict
import pysam
from pysam import Samfile
from collections import Counter
import numpy as np
## class definitions
 
def catch(func, handle=lambda e : e, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return 0

def read_in_vcf(vcf):
  ## function to read in VCF file vcf or files in directory named vcf
  ## this is used to convert snp positions to IDs when merging SNPs
  v = defaultdict(lambda:defaultdict(lambda:defaultdict(str))) # v[chrom][pos] = snpid
  ## check if vcf is file or directory. 
  if os.path.isdir(vcf): ## if directory, fList is list of files in directory to read in
    fList = [os.path.join(vcf,f) for f in os.listdir(vcf)]
  else: ## otherwise, fList has one file to read in
    fList = [vcf]
  for f in fList: ## read in all files in fList
    if f.split('.')[-1]=="gz": ## if file is gzipped, open using gzip library
      vcf_in = gzip.open(f)
    else:
      vcf_in = open(f) 
    sys.stdout.write("Reading in VCF: " + str(f) + "\n")
    for line in vcf_in: ## for each entry in vcf file, store SNP in dictionary
      if line.startswith('#'): # skip over header lines
        continue
      fields = line.rstrip().split()
      chrom = 'chr'+fields[0]
      pos = str(int(fields[1])-1)
      v[chrom][pos] = fields[2]
    vcf_in.close()
  return v


# class definitions for AS Events

class SE:
  # initialized SE object
  def __init__(self,events,readLength,junctionLength,anchorLength):
    self.events = events
    self.readLength = readLength
    self.junctionLength = junctionLength
    self.anchorLength = anchorLength

  def read_in_events(self): ## read in all SE events in AS event file (should be rMATS output format)
    events_in = open(self.events,'r')
    events = defaultdict(lambda:defaultdict(str))
    firstline = True
    header= []
    for line in events_in:
      if firstline:
        header = line.rstrip().split()
        firstline = False
        continue
      ele = line.rstrip().split()
      result = {}
      for i,col in enumerate(header):
        result[col] = ele[i]

      # get effective lengths
      tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
      uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
      dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord

      tLen=min(tE-tS,self.junctionLength/2); ## target exon length
      uLen=min(uE-uS,self.junctionLength/2); ## upstream exon length
      dLen=min(dE-dS,self.junctionLength/2); ## downstream exon length
      
      I_0= min(tLen,self.readLength)+self.readLength-4*self.anchorLength+2; ## effective inclusion form length for JC
      S_0= self.readLength-2*self.anchorLength+1; ## effective skipping form length for JC
      I_1= self.readLength-2*self.anchorLength+tLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target

      ## store event in dictionary
      events[result['chr']][int(result['upstreamES']),int(result['upstreamEE']),int(result['exonStart_0base']),int(result['exonEnd']),int(result['downstreamES']),int(result['downstreamEE'])] = [I_0,S_0,I_1,S_1] + line.rstrip().split()
    return events

class A3SS:
  def __init__(self,events,readLength,junctionLength,anchorLength):
    self.events = events
    self.readLength = readLength
    self.junctionLength= junctionLength
    self.anchorLength =anchorLength

  def read_in_events(self):
    events_in =open(self.events,'r')
    events =defaultdict(lambda:defaultdict(str))
    firstline =True
    header= []
    for line in events_in:
      if firstline:
        header = line.rstrip().split()
        firstline = False
        continue
      ele = line.rstrip().split()
      result = {}
      for i,col in enumerate(header):
        result[col] = ele[i]

      lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
      sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
      fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

      lLen=min(lE-lS,self.junctionLength/2); ## long exon length
      sLen=min(sE-sS,self.junctionLength/2); ## short exon length
      fLen=min(fE-fS,self.junctionLength/2); ## flanking exon length
      aLen=min(lE-lS-(sE-sS),self.junctionLength/2); ## alternative SS region length

      I_0= min(aLen,self.readLength)+self.readLength-4*self.anchorLength+2; ## effective inclusion form length for JC
      S_0= self.readLength-2*self.anchorLength+1; ## effective skipping form length for JC
      I_1= self.readLength-2*self.anchorLength+aLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target


      if result['strand']=="+":
        events[result['chr']][int(result['flankingEE']),int(result['longExonStart_0base']),int(result['flankingEE']),int(result['shortES']),int(result['flankingES']),int(result['shortEE']),'+'] =  [I_0,S_0,I_1,S_1] + line.rstrip().split()
      else:
        events[result['chr']][int(result['longExonEnd']),int(result['flankingES']),int(result['shortEE']),int(result['flankingES']),int(result['shortES']),int(result['flankingEE']),'-'] =  [I_0,S_0,I_1,S_1] + line.rstrip().split()
    return events

class A5SS:
  def __init__(self,events,readLength,junctionLength,anchorLength):
    self.events = events
    self.readLength = readLength
    self.junctionLength = junctionLength
    self.anchorLength =anchorLength

  def read_in_events(self):
    events_in =open(self.events,'r')
    events =defaultdict(lambda:defaultdict(str))
    firstline =True
    header= []
    for line in events_in:
      if firstline:
        header = line.rstrip().split()
        firstline = False
        continue
      ele = line.rstrip().split()
      result = {}
      for i,col in enumerate(header):
        result[col] = ele[i]

      lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
      sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
      fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

      lLen=min(lE-lS,self.junctionLength/2); ## long exon length
      sLen=min(sE-sS,self.junctionLength/2); ## short exon length
      fLen=min(fE-fS,self.junctionLength/2); ## flanking exon length
      aLen=min(lE-lS-(sE-sS),self.junctionLength/2); ## alternative SS region length
      
      I_0= min(aLen,self.readLength)+self.readLength-4*self.anchorLength+2; ## effective inclusion form length for JC
      S_0= self.readLength-2*self.anchorLength+1; ## effective skipping form length for JC
      I_1= self.readLength-2*self.anchorLength+aLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target


      if result['strand']=="+":
        events[result['chr']][int(result['longExonEnd']),int(result['flankingES']),int(result['shortEE']),int(result['flankingES']),int(result['shortES']),int(result['flankingEE']),'+'] =  [I_0,S_0,I_1,S_1] + line.rstrip().split()
      else:
        events[result['chr']][int(result['flankingEE']),int(result['longExonStart_0base']),int(result['flankingEE']),int(result['shortES']),int(result['flankingES']),int(result['shortEE']),'-'] = [I_0,S_0,I_1,S_1] + line.rstrip().split()
    return events

class MXE:
  def __init__(self,events,readLength,junctionLength,anchorLength):
    self.events = events
    self.readLength = readLength
    self.junctionLength= junctionLength
    self.anchorLength =anchorLength

  def read_in_events(self):
    events_in =open(self.events,'r')
    events =defaultdict(lambda:defaultdict(str))
    firstline =True
    header= []
    for line in events_in:
      if firstline:
        header = line.rstrip().split()
        firstline = False
        continue
      ele = line.rstrip().split()
      result = {}

      for i,col in enumerate(header):
        result[col] = ele[i]

      tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
      sS = int(ele[7]); sE = int(ele[8]); ## second exon coord
      uS = int(ele[9]); uE = int(ele[10]); ## upstream exon coord (samller coord)
      dS = int(ele[11]); dE = int(ele[12]); ## downstream exon coord (bigger coord)
  
      tLen=min(tE-tS,self.junctionLength/2); ## target exon length
      sLen=min(sE-sS,self.junctionLength/2); ## second exon length
      uLen=min(uE-uS,self.junctionLength/2); ## upstream exon length
      dLen=min(dE-dS,self.junctionLength/2); ## downstream exon length

      I_0= min(tLen,self.readLength)+self.readLength-4*self.anchorLength+2; ## effective inclusion form length for JC
      S_0= min(sLen,self.readLength)+self.readLength-4*self.anchorLength+2; ## effective skipping form length for JC
      I_1= self.readLength-2*self.anchorLength+tLen+1; ## effective inclusion form length for JC+reads on target
      S_1= self.readLength-2*self.anchorLength+sLen+1; ## effective inclusion form length for JC+reads on target

      events[result['chr']][int(result['upstreamES']),int(result['upstreamEE']),int(result['1stExonStart_0base']),int(result['1stExonEnd']),int(result['2ndExonStart_0base']),int(result['2ndExonEnd']),int(result['downstreamES']),int(result['downstreamEE'])] = [I_0,I_1,S_0,S_1] + line.rstrip().split()
    return events

class RI:
  def __init__(self,events,readLength,junctionLength,anchorLength):
    self.events = events
    self.readLength = readLength
    self.junctionLength = junctionLength
    self.anchorLength = anchorLength

  def read_in_events(self):
    events_in =open(self.events,'r')
    events =defaultdict(lambda:defaultdict(str))
    firstline =True
    header= []
    for line in events_in:
      if firstline:
        header = line.rstrip().split()
        firstline = False
        continue
      ele = line.rstrip().split()
      result = {}
      for i,col in enumerate(header):
        result[col] = ele[i]

      rS = int(ele[5]); rE = int(ele[6]); ## ri exon coord (including up- and down-stream exons)
      uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
      dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord

      rLen=min(rE-rS,self.junctionLength/2); ## ri exon length
      uLen=min(uE-uS,self.junctionLength/2); ## upstream exon length
      dLen=min(dE-dS,self.junctionLength/2); ## downstream exon length
      riLen=min(rE-rS-(uE-uS)-(dE-dS),self.junctionLength/2); ## retained exon length

      I_0= min(riLen,self.readLength)+self.readLength-4*self.anchorLength+2; ## effective inclusion form length for JC
      S_0= self.readLength-2*self.anchorLength+1; ## effective skipping form length for JC
      I_1= self.readLength-2*self.anchorLength+riLen+1; ## effective inclusion form length for JC+reads on target
      S_1= S_0; ## effective skipping form length for JC+reads on target


      events[result['chr']][int(result['upstreamES']),int(result['upstreamEE']),int(result['downstreamES']),int(result['downstreamEE'])] = [I_0,S_0,I_1,S_1] + line.rstrip().split()
    return events

## done definining classes for AS Events

## Class definition for ASAS event
class ASASEvents:
  # initialize object, need list of AS Events for all 5 types, haplotype specific bam files (from rPGA output), and readlengths 
  def __init__(self,a3ss,a5ss,mxe,ri,se,hap1bam,hap2bam,outdir,readLength,junctionLength,anchorLength):
    self.rL = readLength
    self.jL = junctionLength
    self.aL =anchorLength
    self.A3SS = A3SS(a3ss,readLength,junctionLength,anchorLength).read_in_events()
    self.A5SS = A5SS(a5ss,readLength,junctionLength,anchorLength).read_in_events()
    self.MXE = MXE(mxe,readLength,junctionLength,anchorLength).read_in_events()
    self.RI = RI(ri,readLength,junctionLength,anchorLength).read_in_events()
    self.SE = SE(se,readLength,junctionLength,anchorLength).read_in_events()
    self.hap1bam = hap1bam
    self.hap2bam = hap2bam
    self.outdir = outdir


  def get_tags(self,r): ## get read tags 
    return {key: value for key, value in r.tags}

  def is_junction(self,r):
    # check whether r is a junction read     
    for (cigarType,cigarLength) in r.cigar:
      if (cigarType==3):
        return True
    return False

  def get_junction_coordinates(self,r):
    # get junction coordinates that r spans                                                         
    j = [] #list of junctions read spans    
    start_pos = r.pos # read start pos
    genopos = r.pos 
    for (cigarType,cigarLength) in r.cigar:
      if cigarType==0: # matched bases         
        genopos += cigarLength
      elif cigarType==1:# insertion 
        continue
      elif cigarType==2: #deletion  
        genopos += cigarLength
      elif cigarType==3: #skipped bases (junction)
        jStart = genopos #junc start
        jEnd = genopos + cigarLength # junction end
        j.append([int(jStart),int(jEnd)]) # append junction to list j
        genopos += cigarLength
      elif cigarType==4: #soft clipping    
        genopos += cigarLength
      elif cigarType==5: #hard clipping 
        genopos += cigarLength
      else:
        return []
    return j

  def check_mxe_reads(self,readDict,mxeReads,event,isIncl):
    new_reads = []
    for r in mxeReads:
      readID = r[0]
      p = True
      for read in readDict[readID]:
        reference_positions = read.get_reference_positions() # list of positions the read covers 
        if isIncl: # read positiions should just cover exon coordinates agreeing with inclusion isoform
          bad_positions = range(event[1]+1,event[2])+range(event[3]+1,event[6])
        else: # read positions should just cover exon coordinates agreeing wiht exclusion isoform
          bad_positions= range(event[1]+1,event[4])+range(event[4]+1,event[6])
        if any([p for p in reference_positions if p in bad_positions]): # read covers coordinates not agreeing with isoform   
          p = False
          continue
        if self.is_junction(read): # if read is a junctions read, check junctions to make sure they agree with isoform
          jcoords = self.get_junction_coordinates(read) 
          if isIncl: #check junction starts/ends match inclusion form
            for j in jcoords:
              if event[0]<j[0]<event[7] or event[0]<j[1]<event[7]: # jstart or jend falls into isoform region
                if not((j[0]==event[1] and j[1]==event[2]) or (j[0]==event[3] and j[1]==event[6])): # junction doesnt match either us or ds inclusion junction
                  p = False
                  continue
          else: #check junction starts/ends match exclusion form
            for j in jcoords:
              if event[0]<j[0]<event[7] or event[0]<j[1]<event[7]: # jstart or jend falls into isoform region
                if not((j[0]==event[1] and j[1]==event[4]) or (j[0]==event[5] and j[1]==event[6])): # junction doesn't match either us or ds skipping junction
                  p = False
                  continue
      if p: # if both reads in the pair pass, add to new_reads
        new_reads.append(r)
    return new_reads
      

  def check_mxe(self,reads,mxe):
    # check to make sure all reads support mxe events actually support mxe event
    # c_mxe[chrom][a][p][g]['I'][0] = [[rID,rstart,haplotype]]
    for chrom in mxe:
      for a in mxe[chrom]:
        for p in mxe[chrom][a]:
          IJ_R = mxe[chrom][a][p]['R']['I'][0]
          IJ_A = mxe[chrom][a][p]['A']['I'][0]
          SJ_R = mxe[chrom][a][p]['R']['S'][0]
          SJ_A = mxe[chrom][a][p]['A']['S'][0]
          new_ijr = self.check_mxe_reads(reads[chrom],IJ_R,a,True)
          mxe[chrom][a][p]['R']['I'][0] = new_ijr
          new_sjr = self.check_mxe_reads(reads[chrom],SJ_R,a,False)
          mxe[chrom][a][p]['R']['S'][0] = new_sjr
          new_ija = self.check_mxe_reads(reads[chrom],IJ_A,a,True)
          mxe[chrom][a][p]['A']['I'][0] = new_ija
          new_sja = self.check_mxe_reads(reads[chrom],SJ_A,a,False)
          mxe[chrom][a][p]['A']['S'][0] = new_sja
    return 

  def read_in_bam(self,bam,c_a3ss,c_a5ss,c_mxe,c_ri,c_se):
    ## read in allele specific bam file, and assign reads to splice junctions
    sys.stdout.write("Getting junction reads from: " + bam + "\n")
    read_dict = defaultdict(lambda:defaultdict(list)) # read_dict[readID] = [read1,read2]
    bam_in = pysam.Samfile(bam,'rb') ## open bam file as pysam object
    for read in bam_in.fetch(): # for each read in the bam file, determine if it supports an AS Event
      tags = self.get_tags(read) # get read tags
      snppos = tags['SP'].split(';') # snp positions read pair covers
      geno = tags['GT'].split(';') # genotype of each snp (ref or alt)
      haplotype = tags['HT'] #haplotype read is assigned to
      duplicates = [k for k,v in Counter(snppos).items() if v>1] ## get duplicate snps -- means read pair overlaps and both mates cover het snp
      snps = {snppos[i]:geno[i] for i in range(len(snppos))} # dictionary keys=snpposition, value=geno(R or A)
      reference_positions = read.get_reference_positions() # list of positions the read covers
      rstart = mc = read.pos+1 # read start, one based
      rend = reference_positions[-1] # read end position
      chrom = bam_in.getrname(read.reference_id) # chromosome
      rID = read.qname # read id
      read_dict[chrom][rID].append(read)
      if self.is_junction(read): # if read is junction read, check if it supports AS event
        jcoords = self.get_junction_coordinates(read) # junction coordinates of the junction(s) the read covers
        for j in jcoords: # for each junction the read covers
          ## check for SE
          incl = [a for a in self.SE[chrom] if ( 
            (rstart>=a[0] and rend<=a[5]) and 
            ((j[0]==a[1] and j[1]==a[2]) or (j[0]==a[3] and j[1]==a[4])))] ## inclusion form
          skip = [a for a in self.SE[chrom]  if (
            (rstart>=a[0] and rend<=a[5]) and (j[0]==a[1] and j[1]==a[4]))] ## skip form
          for a in incl: # for each AS event the read supports the inclusion form 
            for p,g in snps.iteritems(): # for each SNP, add the read to the count dictionary
              # make sure SNP is not in middle exon
              if ((int(p)<a[2]) or (int(p)>a[3])):
                c_se[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype]) # junction reads only
                c_se[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype]) # junction reads + reads on target
          for a in skip:
            for p,g in snps.iteritems():
              if ((int(p)<a[2]) or (int(p)>a[3])):
                c_se[chrom][a][p][g]['S'][0].append([rID,rstart,haplotype])
                c_se[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
          ## check for MXE
          exon1 = [a for a in self.MXE[chrom]  if (
            (rstart>=a[0] and rend<=a[7]) and
            ((j[0]==a[1] and j[1]==a[2]) or (j[0]==a[3] and j[1]==a[6])))]
          exon1_up = [a for a in self.MXE[chrom]  if ((rstart>=a[0] and rend<=a[7]) and(j[0]==a[1] and j[1]==a[2]))]
          exon1_down = [a for a in self.MXE[chrom]  if ((rstart>=a[0] and rend<=a[7]) and(j[0]==a[3] and j[1]==a[6]))]
          exon2 = [a for a in self.MXE[chrom]  if (
            (rstart>=a[0] and rend<=a[7]) and
            ((j[0]==a[1] and j[1]==a[4]) or (j[0]==a[5] and j[1]==a[6])))]
          exon2_up = [a for a in self.MXE[chrom]  if ((rstart>=a[0] and rend<=a[7]) and(j[0]==a[1] and j[1]==a[4]))]
          exon2_down = [a for a in self.MXE[chrom]  if ((rstart>=a[0] and rend<=a[7]) and(j[0]==a[5] and j[1]==a[6]))]
          for a in exon1:
            for p,g in snps.iteritems():
              if ((int(p)<a[2]) or (int(p)>a[5])):
                c_mxe[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
                c_mxe[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
                if a in exon1_up:
                  c_mxe[chrom][a][p][g]['IU'][0].append([rID,rstart,haplotype])
                if a in exon1_down:
                  c_mxe[chrom][a][p][g]['ID'][0].append([rID,rstart,haplotype])
          for a in exon2:
            for p,g in snps.iteritems():
              if ((int(p)<a[2]) or (int(p)>a[5])):
                c_mxe[chrom][a][p][g]['S'][0].append([rID,rstart,haplotype])
                c_mxe[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
                if a in exon2_up:
                  c_mxe[chrom][a][p][g]['SU'][0].append([rID,rstart,haplotype])
                if a in exon2_down:
                  c_mxe[chrom][a][p][g]['SD'][0].append([rID,rstart,haplotype])
          ## check for ri skipping exon
          skip = [a for a in self.RI[chrom]  if (
            (rstart>=a[0] and rend<=a[3]) and
            (j[0]==a[1] and j[1]==a[2]))]
          for a in skip:
            for p,g in snps.iteritems():
              if ((int(p)<a[1]) or (int(p)>a[2])):
                c_ri[chrom][a][p][g]['S'][0].append([rID,rstart,haplotype])
                c_ri[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
          ## check for A5SS
          longE = [a for a in self.A5SS[chrom]  if (
            (rstart>=a[4] and rend<=a[5]) and
            (j[0]==a[0] and j[1]==a[1]))]
          shortE = [a for a in self.A5SS[chrom]  if (
            (rstart>=a[4] and rend<=a[5]) and
            (j[0]==a[2] and j[1]==a[3]))]
          for a in longE:
            for p,g in snps.iteritems():
              if a[-1]=='+': # pos strand, check SNP doesn't fall in included exon portion
                if ((int(p)<a[2]) or (int(p)>a[0])):
                  c_a5ss[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
                  c_a5ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
              else: # neg strand
                 if ((int(p)<a[1]) or (int(p)>a[3])):
                   c_a5ss[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
                   c_a5ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
          for a in shortE:
            for p,g in snps.iteritems():
              if a[-1]=='+': # pos strand, check SNP doesn't fall in included exon portion
                if ((int(p)<a[2]) or (int(p)>a[0])):
                  c_a5ss[chrom][a][p][g]['S'][0].append([rID,rstart,haplotype])
                  c_a5ss[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
              else: # neg strand
                if ((int(p)<a[1]) or (int(p)>a[3])):
                   c_a5ss[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
                   c_a5ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
          ## check for a3ss
          longE = [a for a in self.A3SS[chrom]  if (
            (rstart>=a[4] and rend<=a[5]) and 
            (j[0]==a[0] and j[1]==a[1]))]
          shortE = [a for a in self.A3SS[chrom]  if (
            (rstart>=a[4] and rend<=a[5]) and
            (j[0]==a[2] and j[1]==a[3]))]
          for a in longE:
            for p,g in snps.iteritems():
              if a[-1]=="+": #pos strand, check SNP doesn't fall in included exon portion  
                if ((int(p)<a[1]) or (int(p)>a[3])):
                  c_a3ss[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
                  c_a3ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
              else: # neg strand
                if ((int(p)<a[2]) or (int(p)>a[0])):
                  c_a3ss[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
                  c_a3ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
          for a in shortE:
            for p,g in snps.iteritems():
              if a[-1]=="+": #pos strand, check SNP doesn't fall in included exon portion
                if ((int(p)<a[1]) or (int(p)>a[3])):
                  c_a3ss[chrom][a][p][g]['S'][0].append([rID,rstart,haplotype])
                  c_a3ss[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
              else:
                if ((int(p)<a[2]) or (int(p)>a[0])):
                  c_a3ss[chrom][a][p][g]['S'][0].append([rID,rstart,haplotype])
                  c_a3ss[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
      else: # if exonic read, check if read supports exon body
        # check if read supports retained intron junction
        incl = [a for a in self.RI[chrom] if (
          (a[0]<=rstart<=a[1] and rend>a[1]) or 
          (rstart<a[2] and a[2]<rend<a[3]))]
        for a in incl:
          for p,g in snps.iteritems():
            if (int(p)<a[1]) or (int(p)>a[2]):
              c_ri[chrom][a][p][g]['I'][0].append([rID,rstart,haplotype])
              c_ri[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
        # check for SE exon reads
        incl = [a for a in self.SE[chrom] if  (a[2]<=rstart and a[3]>=rend)]        
        for a in incl:
          for p,g in snps.iteritems():
            if ((int(p)<a[2]) or (int(p)>a[3])):
              c_se[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
        # check for A3SS reads 
        incl = [a for a in self.A3SS[chrom] if (
            (a[1]<=rstart<=a[3]) and rend<=a[5])]
        for a in incl:
          for p,g in snps.iteritems():
            if a[-1]=="+": #pos strand, check SNP doesn't fall in included exon portion                                                                                                                   
              if ((int(p)<a[1]) or (int(p)>a[3])):
                c_a3ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
            else: #neg strand
              if ((int(p)<a[2]) or (int(p)>a[0])):
                c_a3ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
        # check for A5SS reads
        incl = [a for a in self.A5SS[chrom] if (
            (rstart>a[4]) and (a[2]<=rstart<=a[0]))]
        for a in incl:
          for p,g in snps.iteritems():
            if a[-1]=="+":
              if ((int(p)<a[2]) or (int(p)>a[0])):
                c_a5ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
            else:
              if ((int(p)<a[1]) or (int(p)>a[3])):
                c_a5ss[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
        # check for MXE reads
        incl = [a for a in self.MXE[chrom] if (
            (rstart>=a[2]) and (rend<=a[3]))]
        skip = [a for a in self.MXE[chrom] if (
            (rstart>=a[4]) and (rend<=a[4]))]
        for a in incl:
          for p,g in snps.iteritems():
            if (int(p)<a[2] or int(p)>a[5]):
              c_mxe[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
        for a in skip:
          for p,g in snps.iteritems():
            if (int(p)<a[2] or int(p)>a[5]):
              c_mxe[chrom][a][p][g]['S'][1].append([rID,rstart,haplotype])
        # check for RI reads
        incl = [a for a in self.RI[chrom] if (
            (rstart>=a[1] and rend<=a[2]))]
        for a in incl:
          for p,g in snps.iteritems():
            if (int(p)<a[1] or int(p)>a[2]):
              c_ri[chrom][a][p][g]['I'][1].append([rID,rstart,haplotype])
    bam_in.close()
    self.check_mxe(read_dict,c_mxe)
    return 


  def print_output(self,d,info,outdir,ASType,inclLen,skipLen):
    out1pair = open(os.path.join(outdir,"ASAS.SNP."+ASType+".JunctionReadsOnly.byPair.txt"),'w')
    out1hap = open(os.path.join(outdir,"ASAS.haplotype."+ASType+".JunctionReadsOnly.byPair.txt"),'w')
    out1pair.write('exonID\tREF_IJ\tREF_SJ\tALT_IJ\tALT_SJ\tincLen\tskpLen\n')
    out1hap.write('exonID\tHAP1_IJ\tHAP1_SJ\tHAP2_IJ\tHAP2_SJ\tincLen\tskpLen\n')
    if ASType=="MXE":
      outMXE = open(os.path.join(outdir,"ASAS.SNP.MXE.JunctionReadsOnly.detailed.txt"),"w")
      outMXE.write("exonID\tIJ_REF_UP\tIJ_REF_DOWN\tSJ_REF_UP\tSJ_REF_DOWN\tIJ_ALT_UP\tIJ_ALT_DOWN\tSJ_ALT_UP\tSJ_ALT_DOWN\n")
    for c in d:
      for a in d[c]: # for each as event
        try:
          exonID = info[c][a][4]
        except:
          print c, a
          continue
        IJ = []
        SJ = []
        for p in d[c][a]:
          ID = exonID + '_' + str(c) + '_' + str(p)

          # read counts for read pair as one read
          IJ_Ref2 = len(set(tuple(i)[0] for i in d[c][a][p]['R']['I'][0])) 
          IJ_Alt2 = len(set(tuple(i)[0] for i in d[c][a][p]['A']['I'][0]))
          SJ_Ref2 = len(set(tuple(i)[0] for i in d[c][a][p]['R']['S'][0]))
          SJ_Alt2 = len(set(tuple(i)[0] for i in d[c][a][p]['A']['S'][0]))
            
          out1pair.write(ID + '\t' + str(IJ_Ref2) +'\t' + str(SJ_Ref2) + '\t' + str(IJ_Alt2) + '\t' + str(SJ_Alt2)+'\t' + str(inclLen) + '\t' + str(skipLen) +'\n')

          if ASType=="MXE":
            IJ_RU = len(set(tuple(i)[0] for i in d[c][a][p]['R']['IU'][0]))
            IJ_RD = len(set(tuple(i)[0] for i in d[c][a][p]['R']['ID'][0]))
            SJ_RU = len(set(tuple(i)[0] for i in d[c][a][p]['R']['SU'][0]))
            SJ_RD = len(set(tuple(i)[0] for i in d[c][a][p]['R']['SD'][0]))
            IJ_AU = len(set(tuple(i)[0] for i in d[c][a][p]['A']['IU'][0]))
            IJ_AD = len(set(tuple(i)[0] for i in d[c][a][p]['A']['ID'][0]))
            SJ_AU = len(set(tuple(i)[0] for i in d[c][a][p]['A']['SU'][0]))
            SJ_AD = len(set(tuple(i)[0] for i in d[c][a][p]['A']['SD'][0]))

            outMXE.write(ID + '\t' + str(IJ_RU) +'\t' + str(IJ_RD) + '\t' + str(SJ_RU) + '\t' + str(SJ_RD)+ '\t' + str(IJ_AU) +'\t' + str(IJ_AD) + '\t' + str(SJ_AU) + '\t' + str(SJ_AD)+'\n')

          IJ += d[c][a][p]['R']['I'][0]
          IJ += d[c][a][p]['A']['I'][0]
          SJ += d[c][a][p]['R']['S'][0]
          SJ += d[c][a][p]['A']['S'][0]

        ## accumulate haplotype counts
        IJ_HAP1 = len(set(tuple(i)[0] for i in IJ if i[2]==1))
        SJ_HAP1 = len(set(tuple(i)[0] for i in SJ if i[2]==1))
        IJ_HAP2 = len(set(tuple(i)[0] for i in IJ if i[2]==2))
        SJ_HAP2 = len(set(tuple(i)[0] for i in SJ if i[2]==2))
        ID = exonID
        out1hap.write(ID + '\t' + str(IJ_HAP1) +'\t' + str(SJ_HAP1) + '\t' + str(IJ_HAP2) + '\t' + str(SJ_HAP2)+'\t' + str(inclLen) + '\t' + str(skipLen) +'\n')
        
    out1hap.close()
    out1pair.close()
    if ASType=="MXE":
      outMXE.close()
    return

  def count_splicing_events(self):
    count_se = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(list))))))
    count_ri = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(list))))))
    count_mxe = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(list))))))
    count_a5ss = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(list))))))
    count_a3ss = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(list))))))
    self.read_in_bam(self.hap1bam,count_a3ss,count_a5ss,count_mxe,count_ri,count_se)
    self.read_in_bam(self.hap2bam,count_a3ss,count_a5ss,count_mxe,count_ri,count_se)
    self.print_output(count_se,self.SE,self.outdir,'SE',2,1)
    self.print_output(count_mxe,self.MXE,self.outdir,'MXE',2,2)
    self.print_output(count_ri,self.RI,self.outdir,'RI',2,1)
    self.print_output(count_a3ss,self.A3SS,self.outdir,'A3SS',1,1)
    self.print_output(count_a5ss,self.A5SS,self.outdir,'A5SS',1,1)


class MergeEvents:
  def __init__(self,samples,outdir,vcf,convert,haplotype):
    self.samples = samples
    self.outdir = outdir
    self.vcf = vcf
    self.convert = convert
    self.haplotype = haplotype

  def read_in_samples(self):
    s = []
    samples_in = open(self.samples,'r')
    for line in samples_in:
      s.append(line.rstrip())
    return s

  def add_sample(self,d,l,fn):
    fin = open(fn,'r')
    firstline = True
    for line in fin:
      if firstline:
        firstline = False
        continue
      fields = line.rstrip().split()
      d[fields[0]]['IJC_REF'].append(fields[1])
      d[fields[0]]['SJC_REF'].append(fields[2])
      d[fields[0]]['IJC_ALT'].append(fields[3])
      d[fields[0]]['SJC_ALT'].append(fields[4])
      l[fields[0]] = [fields[5],fields[6]]
    fin.close()
    return

  def add_sample_mxe_detail(self,d,l,fn):
    fin = open(fn,'r')
    firstline = True
    for line in fin:
      if firstline:
        firstline = False
        continue
      fields = line.rstrip().split()
      d[fields[0]]['IJC_REF_U'].append(fields[1])
      d[fields[0]]['IJC_REF_D'].append(fields[2])
      d[fields[0]]['SJC_REF_U'].append(fields[3])
      d[fields[0]]['SJC_REF_D'].append(fields[4])
      d[fields[0]]['IJC_ALT_U'].append(fields[5])
      d[fields[0]]['IJC_ALT_D'].append(fields[6])
      d[fields[0]]['SJC_ALT_U'].append(fields[7])
      d[fields[0]]['SJC_ALT_D'].append(fields[8])
    fin.close()
    return 

  def print_output_mxe_detail_convert(self,d,l,out,v):
     fout = open(out,'w')
     for i in sorted(d):
       event,chrom,pos = i.split('_')
       fout.write( event + '_' + v[chrom][pos]+ '\t'+ ','.join(d[i]['IJC_REF_U'])+ '\t'+ ','.join(d[i]['IJC_REF_D'])+ '\t'+ ','.join(d[i]['SJC_REF_U'])+ '\t'+ ','.join(d[i]['SJC_REF_D'])+ '\t'+ ','.join(d[i]['IJC_ALT_U'])+ '\t'+ ','.join(d[i]['IJC_ALT_D'])+ '\t'+ ','.join(d[i]['SJC_ALT_U'])+ '\t'+ ','.join(d[i]['SJC_ALT_D'])+ '\t' + '\t'.join(l[i])+'\n')
     fout.close()

  def print_output_convert_SNP(self,d,l,out_nofilter,out_filter,v):
    fout_nofilter = open(out_nofilter,'w')
    fout_filter = open(out_filter,'w')
    fout_filter.write('ExonID\tIJC_REF\tSJC_REF\tIJC_ALT\tSJC_ALT\tincLen\tskpLen\n')
    fout_nofilter.write('ExonID\tIJC_REF\tSJC_REF\tIJC_ALT\tSJC_ALT\tincLen\tskpLen\n')
    for i in sorted(d):
      event,chrom,pos = i.split('_')
      if len(d[i]['IJC_REF']) != len(d[i]['SJC_REF']) or len(d[i]['IJC_ALT']) != len(d[i]['SJC_ALT']):
        sys.stderr.write('ERROR' +  str(i) + ' doesn\'t have the same number of replicates for REF or ALT inc/skp\n')
        sys.exit()

      fout_nofilter.write( event + '_' + v[chrom][pos]+ '\t'+ ','.join(d[i]['IJC_REF'])+ '\t'+ ','.join(d[i]['SJC_REF'])+ '\t'+ ','.join(d[i]['IJC_ALT'])+ '\t'+ ','.join(d[i]['SJC_ALT'])+ '\t' + '\t'.join(l[i])+'\n')

      # average psi values for ref and alt
      psi1 = np.average(np.array([catch(lambda: (float(ij)/2)/((float(ij)/2)+float(s))) for ij,s in zip(d[i]['IJC_REF'],d[i]['SJC_REF'])]))
      psi2 = np.average(np.array([catch(lambda: (float(ij)/2)/((float(ij)/2)+float(s))) for ij,s in zip(d[i]['IJC_ALT'],d[i]['SJC_ALT'])]))
      if (min(psi1,psi2)>0.95) or (max(psi1,psi2)<0.05): # psi filters
        continue
      #average readcounts for ref and alt
      rc1 = np.average(np.array([int(ij)+int(s) for ij,s in zip(d[i]['IJC_REF'],d[i]['SJC_REF'])]))
      rc2 = np.average(np.array([int(ij)+int(s) for ij,s in zip(d[i]['IJC_ALT'],d[i]['SJC_ALT'])]))
      if (rc1<10 or rc2<10): # require average read count to be at least 10 for both alleles
        continue
      fout_filter.write(event + '_' + v[chrom][pos]+ '\t'+ ','.join(d[i]['IJC_REF'])+ '\t'+ ','.join(d[i]['SJC_REF'])+ '\t'+ ','.join(d[i]['IJC_ALT'])+ '\t'+ ','.join(d[i]['SJC_ALT'])+ '\t' + '\t'.join(l[i])+'\n')
      
    fout_nofilter.close()
    fout_filter.close()
    return

  def print_output(self,d,l,out_nofilter,out_filter):
    fout_nofilter = open(out_nofilter,'w')
    fout_nofilter.write('ExonID\tIJC_REF\tSJC_REF\tIJC_ALT\tSJC_ALT\tincLen\tskpLen\n')
    fout_filter = open(out_filter,'w')
    fout_filter.write('ExonID\tIJC_REF\tSJC_REF\tIJC_ALT\tSJC_ALT\tincLen\tskpLen\n')

    for i in sorted(d):
      if len(d[i]['IJC_REF']) != len(d[i]['SJC_REF']) or len(d[i]['IJC_ALT']) != len(d[i]['SJC_ALT']):
        sys.stderr.write('ERROR' +  str(i) + ' doesn\'t have the same number of replicates for REF or ALT inc/skp\n')
        sys.exit()

      fout_nofilter.write( str(i)+ '\t'+ ','.join(d[i]['IJC_REF'])+ '\t'+ ','.join(d[i]['SJC_REF'])+ '\t'+ ','.join(d[i]['IJC_ALT'])+ '\t'+ ','.join(d[i]['SJC_ALT'])+ '\t' + '\t'.join(l[i])+'\n')
      # average psi values for ref and alt
      rc1 = np.average(np.array([int(ij)+int(s) for ij,s in zip(d[i]['IJC_REF'],d[i]['SJC_REF'])]))
      rc2 = np.average(np.array([int(ij)+int(s) for ij,s in zip(d[i]['IJC_ALT'],d[i]['SJC_ALT'])]))
      if (rc1<10 or rc2<10): # require average read count to be at least 10 for both alleles
        continue
      psi1 = np.average(np.array([catch(lambda: (float(ij)/2)/((float(ij)/2)+float(s))) for ij,s in zip(d[i]['IJC_REF'],d[i]['SJC_REF'])]))
      psi2 = np.average(np.array([catch(lambda: (float(ij)/2)/((float(ij)/2)+float(s))) for ij,s in zip(d[i]['IJC_ALT'],d[i]['SJC_ALT'])]))
      if (min(psi1,psi2)>0.95) or (max(psi1,psi2)<0.05): # psi filters
        continue
      fout_filter.write( str(i)+ '\t'+ ','.join(d[i]['IJC_REF'])+ '\t'+ ','.join(d[i]['SJC_REF'])+ '\t'+ ','.join(d[i]['IJC_ALT'])+ '\t'+ ','.join(d[i]['SJC_ALT'])+ '\t' + '\t'.join(l[i])+'\n')

    fout_nofilter.close()
    fout_filter.close()
    return 

  def merge_splicing_events(self,ASType):
    samples = self.read_in_samples()
    eventDict1 = defaultdict(lambda:defaultdict(list))
    lengthDict1 = defaultdict(list)
    
    if self.haplotype:
      for s in samples:
        self.add_sample(eventDict1,lengthDict1,os.path.join(s,'ASAS.haplotype.'+ASType+'.JunctionReadsOnly.byPair.txt'))
      self.print_output(eventDict1,lengthDict1,os.path.join(self.outdir,'ASAS.haplotype.'+ASType+'.JunctionReadsOnly.byPair.unfiltered.txt'),os.path.join(self.outdir,'ASAS.haplotype.'+ASType+'.JunctionReadsOnly.byPair.filtered.txt')) 

    else:
      for s in samples:
        self.add_sample(eventDict1,lengthDict1,os.path.join(s,'ASAS.SNP.'+ASType+'.JunctionReadsOnly.byPair.txt'))

      if self.convert:
        self.print_output_convert_SNP(eventDict1,lengthDict1,os.path.join(self.outdir,'ASAS.SNP.'+ASType+'.JunctionReadsOnly.byPair.unfiltered.txt'),os.path.join(self.outdir,'ASAS.SNP.'+ASType+'.JunctionReadsOnly.byPair.filtered.txt'),self.vcf)

      else:
        self.print_output(eventDict1,lengthDict1,os.path.join(self.outdir,'ASAS.SNP.'+ASType+'.JunctionReadsOnly.byPair.unfiltered.txt'),os.path.join(self.outdir,'ASAS.SNP.'+ASType+'.JunctionReadsOnly.byPair.filtered.txt'))

    return

  def merge_mxe_detail(self):
    samples = self.read_in_samples()
    eventDict1 = defaultdict(lambda:defaultdict(list))
    lengthDict1 = defaultdict(list)

    for s in samples:
      self.add_sample_mxe_detail(eventDict1,lengthDict1,os.path.join(s,'ASAS.SNP.MXE.JunctionReadsOnly.detailed.txt'))
    self.print_output_mxe_detail_convert(eventDict1,lengthDict1,os.path.join(self.outdir,'ASAS.SNP.MXE.JunctionReadsOnly.detailed.unfiltered.txt'),self.vcf)
    

def main(args):

  outdir = args.o
  if args.hap1Bam:
    hap1Bam = args.hap1Bam
  else:
    hap1Bam = os.path.join(outdir,'hap1.sorted.bam')
  if args.hap2Bam:
    hap2Bam = args.hap2Bam
  else:
    hap2Bam = os.path.join(outdir,'hap2.sorted.bam')

  if args.asdir:
    se_fn = os.path.join(args.asdir,"fromGTF.SE.txt")
    mxe_fn = os.path.join(args.asdir,"fromGTF.MXE.txt")
    ri_fn = os.path.join(args.asdir,"fromGTF.RI.txt")
    a5ss_fn = os.path.join(args.asdir,"fromGTF.A5SS.txt")
    a3ss_fn = os.path.join(args.asdir,"fromGTF.A3SS.txt")
  else:
    se_fn = args.SE
    mxe_fn = args.MXE
    ri_fn = args.RI
    a3ss_fn = args.A3SS
    a5ss_fn = args.A5SS

  samples_fn = args.samples

  if args.merge:
    if not args.samples:
      sys.stderr.write('rPGA: ERROR! Need to provide file containing list of sample output directories when using rPGA splicing --merge\n\n')
      sys.exit()
    if not os.path.exists(outdir):
      os.mkdir(outdir)
    
    if args.mxeDetail:
      if not args.v:
        sys.stderr.write('rPGA: ERROR! Need to provide a vcf file to convert snp positions to snp IDs\n')
        sys.exit()
      if not os.path.exists(outdir):
        os.mkdir(outdir)
      convert = True
      VCF = read_in_vcf(args.v)
      p = MergeEvents(samples_fn,outdir,VCF,convert,False)
      p.merge_mxe_detail()

    elif args.pos2id:
      if not args.v:
        sys.stderr.write('rPGA: ERROR! Need to provide a vcf file to convert snp positions to snp IDs\n')
        sys.exit()
      if not os.path.exists(outdir):
        os.mkdir(outdir)
      convert = True
      VCF = read_in_vcf(args.v)
      p = MergeEvents(samples_fn,outdir,VCF,convert,args.haplotype)
      p.merge_splicing_events('SE')
      p.merge_splicing_events('A3SS')
      p.merge_splicing_events('A5SS')
      p.merge_splicing_events('MXE')
      p.merge_splicing_events('RI')

    else:
      convert = False
      p = MergeEvents(samples_fn,outdir,args.v,convert,args.haplotype)
      p.merge_splicing_events('SE')
      p.merge_splicing_events('A3SS')
      p.merge_splicing_events('A5SS')
      p.merge_splicing_events('MXE')
      p.merge_splicing_events('RI')

  else:
    if args.readlength:
      readLength = int(args.readlength)
    else:
      readLength = 100
    if args.anchorlength:
      anchorLength = int(args.anchorlength)
    else:
      anchorLength = 8
    junctionLength = 2*(int(readLength) - int(anchorLength))
    p = ASASEvents(a3ss_fn,a5ss_fn,mxe_fn,ri_fn,se_fn,hap1Bam,hap2Bam,outdir,readLength,junctionLength,anchorLength)
    p.count_splicing_events()
    
