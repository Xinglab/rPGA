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
import sys, os
import subprocess
import re,logging,time,datetime,commands,argparse
from collections import defaultdict
import pysam, pybedtools
from pysam import Samfile
from pybedtools import BedTool

################################################################################
##            PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                ##
################################################################################

def STAR_create_genome(project, genome, gnme, threads):
  # build up options
  opts = ""
  opts += (" --runMode genomeGenerate")  
  opts += (" --genomeDir " + str(project) + "/" + str(gnme) + "/" + "STARindex")
  opts += (" --genomeFastaFiles " + str(genome))
  opts += (" --runThreadN "+str(threads)) 

  env_cpy = os.environ.copy()
  commandSTAR = ("STAR" + " " + opts)

  oFile = open(str(project) + "/mapping_commands.sh","w")
  oFile.write("##### Creating Genome for " + str(gnme) + "#####\n" +\
              commandSTAR + "\n#\n")
  oFile.flush()
  oFile.close()
  status,output=commands.getstatusoutput(commandSTAR)

  return


def STAR_perform_mapping(project, gnme, seqs, threads,mismatches,gz, multimapped):

  # build up options
  opts = ""
  opts += (" --genomeDir " + str(project) + "/" + str(gnme) + "/" + "STARindex --readFilesIn ")
  for line in open(seqs) :
    line = line.rstrip()
    opts += (str(line) + " ")
  if gz:
    opts += (" --readFilesCommand gunzip -c")
  opts += (" --runThreadN " + str(threads))
  opts += (" --outFilterMultimapNmax " + str(multimapped))
  opts += (" --alignEndsType EndToEnd")
  opts += (" --outFilterMismatchNmax " + str(mismatches))
  opts += (" --outFileNamePrefix " + str(project))
  opts += ("/" + str(gnme) + "/" + "STARalign/ ")
  opts += ("--outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonical")
  opts += (" --alignIntronMax 300000 --outSJfilterOverhangMin -1 8 8 8")



  env_cpy = os.environ.copy()
  commandSTAR = ("STAR" + " " + opts)
  print commandSTAR
  oFile = open(str(project) + "/mapping_commands.sh","w")
  oFile.write("##### Creating Genome for " + str(gnme) + "#####\n" +\
              commandSTAR + "\n#\n")
  oFile.flush()
  oFile.close()
  status,output=commands.getstatusoutput(commandSTAR)

  return


def sam_to_sorted_bam(sam_fn):
  print "sam to bam: ", sam_fn
  bam_fn = sam_fn + '.bam'
  pysam.view("-Sb","-o%s" % bam_fn, sam_fn+".sam")
  pysam.sort(bam_fn, sam_fn+'.sorted')
  pysam.index(sam_fn+'.sorted.bam')
  return

def worker(i):
  p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, i, writeBam,discoverJunctions,rnaedit,editFile)
  p.haplotype_specific_junctions()

################################################################################
#             PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                 #
################################################################################

class PersonalizeGenome :
  def __init__(self, outDir, vcf, ref, hap1Ref, hap2Ref,rnaedit,editFile):
    self._outDir = outDir
    self._vcf = vcf
    self._ref = ref
    self._hap1Ref = hap1Ref
    self._hap2Ref = hap2Ref
    self._rnaedit = rnaedit
    self._editFile = editFile
    self._report = os.path.join(outDir,'report.personalize.txt')

  def read_reference(self):
    f = defaultdict(list)
    ref_in = open(self._ref)
    key = ''
    for line in ref_in:
      line = line.rstrip()
      if line.startswith('>'):
        key = line[1:]
      else:
        for b in line:
          f[key].append(b)
    ref_in.close()
    return f

  def read_in_edit(self):
    edit_in = open(self._editFile)
    e = defaultdict(list)
    firstline = True
    for line in edit_in:
      if firstline:
        header = line.rstrip().split()
        firstline = False
      else:
        fields = line.rstrip().split()
        result = {}
        for i,col in enumerate(header):
          result[col] = fields[i]
        e[result['chromosome'][3:]].append(int(result['position']))
    edit_in.close()
    return e

  def read_in_vcf(self):
    CHROMS = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    VCF_HEADER = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    v1,v2 = defaultdict(lambda: defaultdict(list)),defaultdict(lambda: defaultdict(list))
    for c in CHROMS:
      vcf_in = open(self._vcf + '/' + c +'.vcf')
      for line in vcf_in:
        if line.startswith('#'):
          continue
        else:
          result = {}
          fields = line.rstrip().split()
          for i,col in enumerate(VCF_HEADER):
            result[col] = fields[i]
          infos = [x for x in result['INFO'].split(';') if x.strip()]
          for i in infos:
            if '=' in i:
              key,value = i.split('=')
              result[key] = value
          if (result['VT'] == 'SNP'): 
            geno = result['SAMPLE'].split(':')[0]
            if '|' in geno:
              g1 = int(geno.split('|')[0])
              g2 = int(geno.split('|')[1])
              alt = result['ALT'].split(',')
              alleles = [result['REF']] + alt
              if ( re.match(r'[ACGT]',alleles[0]) and re.match(r'[ACGT]',alleles[g1]) and re.match(r'[ACGT]',alleles[g2])):
                if (g1 != 0):
                  v1[result['CHROM']][int(result['POS'])] = [alleles[0],alleles[g1]]
                if (g2 != 0):
                  v2[result['CHROM']][int(result['POS'])] = [alleles[0],alleles[g2]]
    vcf_in.close()
    return v1,v2

  def personalize_genome(self):
    if self._rnaedit:
      hap1 = self.read_reference()
      vcf1,vcf2 = self.read_in_vcf()
      editSites = self.read_in_edit()
      report_out = open(self._report,'w')
      hap1_out = open(self._hap1Ref,'w')
      editCounter = 0
      snpCounter = 0
      for chrom in vcf1:
        for pos in vcf1[chrom]:
          if pos in editSites[chrom]:
            editCounter += 1
            hap1['chr'+chrom][int(pos)-1] = 'N'
          else:
            snpCounter += 1
            ref = vcf1[chrom][pos][0]
            alt = vcf1[chrom][pos][1]
            hap1['chr'+chrom][int(pos)-1] = alt
      for chrom in hap1:
        hap1_out.write('>'+chrom+'\n')
        hap1_out.write(''.join(hap1[chrom])+'\n')
      hap1.clear()
      hap1_out.close()

      report_out.write('# number of hap1 SNPs overlapping RNA editing sites: '+ str(editCounter) + '\n')
      report_out.write('# number of hap1 SNPs not overlapping RNA editing sites' + str(snpCounter) + '\n')
      report_out.write('# total number of changed bases in hap1: '+str(editCounter + snpCounter) + '\n')

      editCounter = 0
      snpCounter = 0
      hap2 = self.read_reference()
      hap2_out = open(self._hap2Ref,'w')
      for chrom in vcf2:
        for pos in vcf2[chrom]:
          if pos in editSites[chrom]:
            editCounter+=1
            hap2['chr'+chrom][int(pos)-1] = 'N'
          else:
            snpCounter += 1
            ref = vcf2[chrom][pos][0]
            alt = vcf2[chrom][pos][1]
            hap2['chr'+chrom][int(pos)-1] = alt
      for chrom in hap2:
        hap2_out.write('>'+chrom+'\n')
        hap2_out.write(''.join(hap2[chrom])+'\n')
      hap2.clear()
      hap2_out.close()
      
      report_out.write('# number of hap2 SNPs overlapping RNA editing sites: '+ str(editCounter) + '\n')
      report_out.write('# number of hap2 SNPs not overlapping RNA editing sites' + str(snpCounter) + '\n')
      report_out.write('# total number of changed bases in hap2: '+str(editCounter + snpCounter) + '\n')

    else:
      hap1 = self.read_reference()
      vcf1,vcf2 = self.read_in_vcf()
      report_out = open(self._report,'w')
      hap1_out = open(self._hap1Ref,'w')
      snpCounter = 0
      for chrom in vcf1:
        for pos in vcf1[chrom]:
          snpCounter += 1
          ref = vcf1[chrom][pos][0]
          alt = vcf1[chrom][pos][1]
          hap1['chr'+chrom][int(pos)-1] = alt
      for chrom in hap1:
        hap1_out.write('>'+chrom+'\n')
        hap1_out.write(''.join(hap1[chrom])+'\n')
      hap1.clear()
      hap1_out.close()
      report_out.write('# number of hap1 SNPs: '+ str(snpCounter) + '\n')
      
      snpCounter = 0
      hap2 = self.read_reference()
      hap2_out = open(self._hap2Ref,'w')
      for chrom in vcf2:
        for pos in vcf2[chrom]:
          snpCounter += 1
          ref = vcf2[chrom][pos][0]
          alt = vcf2[chrom][pos][1]
          hap2['chr'+chrom][int(pos)-1] = alt
      for chrom in hap2:
        hap2_out.write('>'+chrom+'\n')
        hap2_out.write(''.join(hap2[chrom])+'\n')
      hap2.clear()
      hap2_out.close()
      report_out.write('# number of hap2 SNPs: '+ str(snpCounter) + '\n')

    return

class DiscoverSpliceJunctions :
  def __init__(self, outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting,rnaedit,editFile):
    self._outDir = outDir
    self._vcf = os.path.join(vcf,chromosome+'.vcf') 
    self._gtf = gtf
    self._hap1Bam = hap1Bam
    self._hap2Bam = hap2Bam
    self._refBam = refBam
    self._chromosome = chromosome
    self._writeBam = writeBam
    self._discoverJunctions = discoverJunctions
    self._conflicting = writeConflicting
    self._VCFheader = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    self._Bam1Out = outDir + '/hap1.' + chromosome + '.as.bam'
    self._Bam2Out = outDir + '/hap2.' + chromosome + '.as.bam'
    self._rnaedit = rnaedit
    self._editFile = editFile
    self._report = os.path.join(outDir, "report."+chromosome+".txt")


  def read_in_rna_editing(self):
    logging.info('reading in rna editing events')
    e = defaultdict(list)
    fin = open(self._editFile)
    header = []
    firstline = True
    for line in fin:
      if firstline:
        header = line.rstrip().split()
        firstline = False
      else:
        fields = line.rstrip().split()
        result = {}
        for i,col in enumerate(header):
          result[col] = fields[i]
          e[result['chromosome'][3:]].append(int(result['position'])-1)
    return e

  def read_in_vcf(self):
    vcf_in = open(self._vcf)
    v = defaultdict(list)
    vids = defaultdict(str)

    if self._rnaedit:
      editPos = self.read_in_rna_editing()
      for line in vcf_in:
        if line.startswith('#'): # skip over header lines
          continue
        else:
          result = {} # store line in dictionary
          fields = line.rstrip().split()
          for i,col in enumerate(self._VCFheader):
            result[col] = fields[i]
          infos = [x for x in result['INFO'].split(';') if x.strip()]
          for i in infos:
            if '=' in i:
              key,value = i.split('=')
              result[key] = value
          if (int(result['POS'])-1) in editPos[result['CHROM']]:
            continue
          else:
            if (result['VT'] == 'SNP'):
              geno = result['SAMPLE'].split(':')[0] 
              if '|' in geno: ## if genotype is phased
                g1 = int(geno.split('|')[0]) ## hap1 genotype
                g2 = int(geno.split('|')[1]) ## hap2 genotype
                alt = result['ALT'].split(',') ## alternate allele(s) in a list
                alleles = [result['REF']] + alt ## [ref, alt1, alt2,...]
                if ( re.match(r'[ACGT]',alleles[0]) and re.match(r'[ACGT]',alleles[g1]) and re.match(r'[ACGT]',alleles[g2])): ## make sure ref and alt alleles are A,C,T,or G
                  if (g1!=g2): # heterozygous snp
                    v[int(result['POS'])-1] = [alleles[g1],alleles[g2]] ## subtract one from position (1 based) to match bam file (0 based)
                    vids[int(result['POS'])-1] = result['ID']
                  else:
                    vids[int(result['POS'])-1] = result['ID']


    else:
      for line in vcf_in:
        if line.startswith('#'): # skip over header lines 
          continue
        else:
          result = {} # store line in dictionary
          fields = line.rstrip().split()
          for i,col in enumerate(self._VCFheader):
            result[col] = fields[i]
          infos = [x for x in result['INFO'].split(';') if x.strip()]
          for i in infos:
            if '=' in i:
              key,value = i.split('=')
              result[key] = value
          if (result['VT'] == 'SNP'):
            geno = result['SAMPLE'].split(':')[0]
            if '|' in geno: ## if genotype is phased                                                                                                                                                      
              g1 = int(geno.split('|')[0]) ## hap1 genotype                                                                                                                                               
              g2 = int(geno.split('|')[1]) ## hap2 genotype                                                                                                                                               
              alt = result['ALT'].split(',') ## alternate allele(s) in a list                                                                                                                             
              alleles = [result['REF']] + alt ## [ref, alt1, alt2,...]                                                                                                                                    
              if ( re.match(r'[ACGT]',alleles[0]) and re.match(r'[ACGT]',alleles[g1]) and re.match(r'[ACGT]',alleles[g2])): ## make sure ref and alt alleles are A,C,T,or G                               
                if (g1!=g2): # heterozygous snp                                                                                                                                                           
                  v[int(result['POS'])-1] = [alleles[g1],alleles[g2]] ## subtract one from position (1 based) to match bam file (0 based)                                                                 
                  vids[int(result['POS'])-1] = result['ID']
                else:
                  vids[int(result['POS'])-1] = result['ID']

    vcf_in.close()
    return v,vids
    
  def read_in_gtf(self):
    gtf_in = open(self._gtf)
    GTF_HEADER  = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
    gg = defaultdict(lambda:defaultdict(list)) #geneGroup[chrom][pos] = [transcript ids]
    g = defaultdict(lambda:defaultdict(list)) #gtf[chrom][transcript]=[list of exons in transcript]
    gi = defaultdict(list) # geneInfo[transcript_id] = [gene_id,gene_name]
    for line in gtf_in:
      if line.startswith('#'):
        continue
      else:
        result = {}
        fields = line.rstrip().split('\t')
      for i,col in enumerate(GTF_HEADER):
        result[col] = fields[i]
      infos = [x.strip() for x in re.split(';',fields[8]) if x.strip()]
      for i in infos:
        key,value = re.split(' ',i)
        result[key] = value
      if result['feature']=="transcript":
        g[result['chrom']][result['transcript_id']] = []
        group = range(int(result['start'])/1000,int(result['end'])/1000 + 1)
        for i in group:
          gg[result['chrom']][i].append(result['transcript_id'])
          gi[result['transcript_id']] = [result['gene_id'],result['gene_name'],result['strand']]
      elif result['feature']=='exon':
        g[result['chrom']][result['transcript_id']].append([int(result['start'])-1,int(result['end'])+1])
    gtf_in.close()
    return gg,g,gi

  def characterize_junction(self,chrom,start,end,gg,g,gi):
    ## characterize junctions as novel or reference
    chrom = 'chr'+chrom
    novel3 = True
    novel5 = True
    NC = True
    s = "N"
    strand = '+'
    transcripts = gg[chrom][start/1000]
    for t in transcripts:
      strand = gi[t][2]
      if len(g[chrom][t])>1: # at least 2 exons in transcript
        g[chrom][t] = sorted(g[chrom][t])
        for i in range(len(g[chrom][t])-1): #check every splice junction in transcript to see if SJ in transcript
          if g[chrom][t][i][1]==start:
            if strand=="+":
              novel5=False
              if g[chrom][t][i+1][0]==end:
                novel3=False
                NC=False
                break
            else:
              novel3=False
              if g[chrom][t][i+1][0]==end:
                novel5=False
                NC=False
                break
          elif g[chrom][t][i+1][0]==end:
            if strand=="+":
              novel3=False
            else:
              novel5=False
            if not NC:
              break
    if (novel3 and  novel5):
      s = "N35"
    elif novel3:
      s = "N3"
    elif novel5:
      s = "N5"
    elif NC:
      s = "NC"
    else:
      s = "R"
    return s,strand 

  def get_junction_coordinates(self,r):
    # get junction coordinates that r spans
    j = [] #list of junctions read spans
    start_pos = r.pos
    genopos = r.pos
    for (cigarType,cigarLength) in r.cigar:
      if cigarType==0: # matched bases
        genopos += cigarLength
      elif cigarType==1:# insertion
        continue
      elif cigarType==2: #deletion
        genopos += cigarLength
      elif cigarType==3: #skipped bases (junction)                            
        jStart = genopos+1
        jEnd = genopos + cigarLength
        j.append([int(jStart),int(jEnd)])
        genopos += cigarLength
      elif cigarType==4: #soft clipping
        genopos += cigarLength
      elif cigarType==5: #hard clipping
        genopos += cigarLength
      else: 
        return []
    return j


  def check_cigar(self,r, allele, p):
    # checks if read r matches allele a at position p 
    # return 0 - r is allele specific 
    # return 1 - snp falls into intron
    # return 2 - read has non N/M in cigar
    # return 3 - read is not allele specific  
    readpos = 0 # keeps track of position in read
    genopos = r.pos #keeps track of position in genome
    for cigar in  r.cigar: # go through each cigar tuple at a time                              
      if cigar[0]==0: #match
      # check if snp in this block of matched bases           
        if p <= ( genopos + int(cigar[1])): 
          readbase = r.seq[p - genopos + readpos] # check the read nt at the position of the snp
          if readbase == allele: # if the read matches allele at p, return 0 (allele specific)
            return 0
          else:
            return 3
        else:
          genopos += int(cigar[1]) # otherwise, continue searching throught the rest of the read
          readpos += int(cigar[1])
      elif cigar[0]==3: # intron                     
        genopos += int(cigar[1])
        if p <= (genopos): # snp falls into intron
          return 1
      else: #read has non N/M in cigar, throw it out                                              
        return 2
    return 3 # read is not allele specific

  def check_variant_in_read(self,r,p):
    readpos = 0
    genopos = r.pos
    for cigar in r.cigar:
      if cigar[0]==0:
        if p <= ( genopos + int(cigar[1])):
          return 0
        else:
          genopos += int(cigar[1])
          readpos += int(cigar[1])
      elif cigar[0]==3:
        genopos += int(cigar[1])
        if p <= (genopos): 
          return 1
      else:
        return 1
    return 1


  def num_mismatches(self,r):
    # return number of mismatches
    nm = 0
    for tag in r.tags:
      if tag[0]=='nM': # nM is tag for mismatches
        return int(tag[1])
    return nm
  

  def is_unique(self,r):
    # check if read is uniquely mapped
    for tag in r.tags:
      if tag[0]=='NH': # NH is tag for number of hits
        if int(tag[1])==1: # uniquely mapped if NH=1
          return True
    return False

  def is_junction(self,r):
    # check whether r is a junction read
    for (cigarType,cigarLength) in r.cigar:
      if (cigarType==3):
        return True
    return False

  def get_splicesite_snp(self,start,end,v):
    # returns splice site snp in junction start-end
    variant_positions = [p for p in (range(start-1,start+1)+range(end-2,end)) if p in v]
    if len(variant_positions)>0:
      return [v[i] for i in variant_positions]
    else:
      return []

  def haplotype_specific_read(self,r,v,h):
    # return 0 = hap spec, 1 = not hap spec, 2 = conflicting, 3 = multimapped
    if not self.is_unique(r):
      return 3
    readpos = 0
    genopos = r.pos
    snps = []
    for (cigarType,cigarLength) in r.cigar:
      if cigarType==0:
        snps += [p for p in range(genopos,genopos+int(cigarLength)) if p in v]
        genopos += int(cigarLength)
        readpos += int(cigarLength)
      elif cigarType==3:
        genopos += cigarLength
      else:
        return 1
    if len(snps)==0:
      return 1
    else:
      check = [self.check_cigar(r, v[p][h], p) for p in snps]
      if check.count(0)>check.count(3):
        return 0
      elif check.count(0)==check.count(3):
        return 2
      else:
        return 1

      

  def haplotype_specific_junctions(self):
    bam1 = pysam.Samfile(self._hap1Bam)
    bam2 = pysam.Samfile(self._hap2Bam)
    snpreads1,snpreads2 = defaultdict(list),defaultdict(list)
    reads1,reads2 = defaultdict(list),defaultdict(list)
    hetsnps,snpids = self.read_in_vcf()
    geneGroup,gtf,geneInfo = self.read_in_gtf()
    spec1,spec2  = list(),list()

    for r in bam1.fetch('chr'+str(self._chromosome)):
      spec = self.haplotype_specific_read(r,hetsnps,0)
      snpreads1[spec].append(r.qname)
      reads1[r.qname].append(r)

    for r in bam2.fetch('chr'+str(self._chromosome)):
      spec = self.haplotype_specific_read(r,hetsnps,1)
      snpreads2[spec].append(r.qname)
      reads2[r.qname].append(r)

    conflicting = list(set([ r for r in snpreads1[0] if r in snpreads2[0]] + snpreads1[2] + snpreads2[2]))
    snpreads1[0] = list(set([ r for r in snpreads1[0] if ((r not in conflicting) and (r not in snpreads2[3]))]))
    snpreads2[0] = list(set([ r for r in snpreads2[0] if ((r not in conflicting) and (r not in snpreads1[3]))]))

    for qname in snpreads1[0]:
      if qname in reads2:
        if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads1[qname]))]): # reads have same starting position in hap1 and hap2
          if self.num_mismatches(reads1[qname][0]) < self.num_mismatches(reads2[qname][0]):
            spec1.append(qname)
    
    for qname in snpreads2[0]:
      if qname in reads1:
        if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads2[qname]))]): # reads have same starting position in hap1 and hap2 
          if self.num_mismatches(reads2[qname][0]) < self.num_mismatches(reads1[qname][0]):
            spec2.append(qname)

    conflictCount = len(set(conflicting))
    hap1Count = len(set(spec1))
    hap2Count = len(set(spec2))
                    
    report_out = open(self._outDir + '/report.'+self._chromosome+'.txt','w')
    report_out.write('# ' + self._chromosome + '\t' + 'hap1' + '\t' + str(hap1Count) + '\n')
    report_out.write('# ' + self._chromosome + '\t' + 'hap2' +'\t' + str(hap2Count) + '\n')
    report_out.write('# ' + self._chromosome + '\t' + 'conflicting' +'\t' + str(conflictCount) + '\n')


    if self._writeBam:
      out1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.bam','wb',template=bam1)
      out2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.bam','wb',template=bam2)
      for qname in spec1:
        for r in reads1[qname]:
          out1.write(r)
      for qname in spec2:
        for r in reads2[qname]:
          out2.write(r)

    if self._conflicting:
      conflict1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.conflicting.bam','wb',template=bam1)
      conflict2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.conflicting.bam','wb',template=bam2)
      for qname in conflicting:
        for r in reads1[qname]:
          conflict1.write(r)
        for r in reads2[qname]:
          conflict2.write(r)
    
    if self._discoverJunctions:
      bamr = pysam.Samfile(self._refBam,"rb")
      junctions = defaultdict(lambda: defaultdict(set))
      for qname in reads1:
        if ((qname not in spec2) and (qname not in conflicting)):
          for r in reads1[qname]:
            juncs = self.get_junction_coordinates(r)
            for j in juncs:
              start,end = j
              junctions['1'][start,end].add(r.pos)

      for qname in reads2:
        if ((qname not in spec1) and (qname not in conflicting)):
          for r in reads2[qname]:
            juncs = self.get_junction_coordinates(r)
            for j in juncs:
              start,end = j
              junctions['2'][start,end].add(r.pos)

      for r in bamr.fetch('chr'+self._chromosome):
        if (r.qname not in conflicting):
          juncs = self.get_junction_coordinates(r)
          for j in juncs:
            start,end = j
            junctions['R'][start,end].add(r.pos)

      nR = {}
      for i in ['1','2','R']:
        nR[i] = {j:len(junctions[i][j]) for j in junctions[i]}
    
      spec = {}
      spec['1'] = [j for j in nR['1'] if ((j not in nR['2']) and (j not in nR['R']) and nR['1'][j]>1 and len(self.get_splicesite_snp(j[0],j[1],snpids))>0)]
      spec['2'] = [j for j in nR['2'] if ((j not in nR['1']) and (j not in nR['R']) and nR['2'][j]>1 and len(self.get_splicesite_snp(j[0],j[1],snpids))>0)]
      spec['12'] = [j for j in nR['1'] if ((j in nR['2']) and (j not in nR['R']) and nR['1'][j]>1 and nR['2'][j]>1 and len(self.get_splicesite_snp(j[0],j[1],snpids))>0)]
      spec['R'] = [j for j in nR['R'] if ((j not in nR['1']) and (j not in nR['2'])and nR['R'][j]>1 and len(self.get_splicesite_snp(j[0],j[1],snpids))>0)]
      bed = defaultdict(list)
      for h in ['1','2','R']:
        counter = 0
        for start,end in spec[h]:
          counter += 1
          num_overlapping_reads = sum([nR[h][j] for j in nR[h] 
                                     if ( j[0]<end and j[1]>start and
                                          self.characterize_junction(self._chromosome,j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
          if num_overlapping_reads > 0:
            freq = float(nR[h][start,end])/float(num_overlapping_reads)
          else:
            freq = 1
          n_or_r,strand = self.characterize_junction(self._chromosome,start,end,geneGroup,gtf,geneInfo)
          snp = ','.join(self.get_splicesite_snp(start,end,snpids))
          bed[h].append('chr'+self._chromosome+' '+str(start) + ' ' + str(end) + ' J_'+str(counter)+'_'+n_or_r+'_'+snp+ ' ' + str(strand) + ' ' +str(freq))
          
      counter = 0
      for start,end in spec['12']:
        counter += 1
        num_overlapping_reads1 = sum([nR['1'][j] for j in nR['1'] if 
                                    ( j[0]<end and j[1]>start and 
                                      self.characterize_junction(self._chromosome,j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
        num_overlapping_reads2 = sum([nR['2'][j] for j in nR['2'] if 
                                    ( j[0]<end and j[1]>start and 
                                      self.characterize_junction(self._chromosome,j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
        if num_overlapping_reads1 == 0:
          freq1=1.0
        else:
          freq1 = float(nR['1'][start,end])/float(num_overlapping_reads1)
        if num_overlapping_reads2==0:
          freq2 = 1.0
        else:
          freq2 = float(nR['2'][start,end])/float(num_overlapping_reads2)
        n_or_r,strand = self.characterize_junction(self._chromosome,start,end,geneGroup,gtf,geneInfo)
        snp =','.join(self.get_splicesite_snp(start,end,snpids))
        bed['12'].append('chr'+self._chromosome+' '+str(start) + ' ' + str(end) + ' J_'+str(counter)+'_'+n_or_r+'_'+snp+ ' ' + str(strand)+' '+str((freq1+freq2)/2))
      files_dict = {}
      files_dict['1'] = self._outDir + '/hap1.'+self._chromosome+'.specific.bed'
      files_dict['2'] = self._outDir + '/hap2.'+self._chromosome+'.specific.bed'
      files_dict['12'] = self._outDir + '/hap1hap2.'+self._chromosome+'.specific.bed'
      files_dict['R'] = self._outDir + '/ref.'+self._chromosome+'.specific.bed'
      for h in ['1','2','12','R']:
        bedstring = '\n'.join(bed[h])
        mybed = BedTool(bedstring, from_string=True).sort().saveas(files_dict[h])
    return


### end of defining functions and classes ##

def isHelpString(s) :
  norm = s.strip().lower()
  return norm == "help" or norm == "--help"

def main(args) :
  """
    Main entry point for this script.

    :param args: the arguments for this script, as a list of string. Should
                 already have had things like the script name stripped. That
                 is, if there are no args provided, this should be the empty
                 list.
  """

  helpStr = "------------------------------------------------------------\n" +\
            "                        RUNNING rPGA                        \n" +\
            "------------------------------------------------------------\n" +\
            "To run personalize, where $ is your prompt:                 \n" +\
            "                                                            \n" +\
            "$ rPGA run personalize                                      \n" +\
            "                                                            \n" +\
            "To run mapping, where $ is your prompt:                     \n" +\
            "                                                            \n" +\
            "$ rPGA run mapping                                          \n" +\
            "                                                            \n" +\
            "If gzipped read sequence file(s):                           \n" +\
            "                                                            \n" +\
            "$ rPGA run mapping -g                                       \n" +\
            "                                                            \n" +\
            "If you are planning to run rPGA alleles and only need       \n" +\
            "alignments to hap1 and hap2 personal genomes, run:          \n" +\
            "                                                            \n" +\
            "$ rPGA run mapping alleles -g                               \n" +\
            "                                                            \n" +\
            "To run discover, where $ is your prompt:                    \n" +\
            "                                                            \n" +\
            "$ rPGA run discover -c CHROM                                \n" +\
            "                                                            \n" +\
            "To run alleles, where $ is your prompt:                     \n" +\
            "                                                            \n" +\
            "$ rPGA run alleles -c CHROM                                 \n" 


  command = args.command
  if args.T:
    threads = int(args.T)
  else:
    threads = 8

  if args.N:
    mismatches = int(args.N)
  else:
    mismatches = 3
     
  if args.M:
    multimapped = int(args.M)
  else:
    multimapped = 20

  if args.rnaedit:
    editFile = args.e
  else:
    editFile = ""


  chromosome = args.c
  writeBam = args.b
  multiprocessing = args.p
  gzipped = args.g
  writeConflicting = args.conflict
  rnaedit = args.rnaedit

  if len(command) == 1 or (len(command)==2 and isHelpString(command[1].strip().lower())):
    sys.stderr.write(helpStr + "\n\n")
  else :
    setting = command[1].strip().lower()
    outDir = open(".rPGAProject.yaml").readline().rstrip()
    hap1Ref = os.path.join(outDir, "hap1.fa")
    hap2Ref = os.path.join(outDir, "hap2.fa")
    vcf = ""
    gtf = ""
    ref = ""
    seqs = ""

    if setting == "personalize" :
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2 :
        sys.stderr.write("Input arguments are not correct\n")
        sys.stderr.write(helpStr + '\n\n')
        sys.exit()
      else :
        ref = open(".rPGAGenome.yaml").readline().rstrip()
        vcf = open(".rPGAGenotype.yaml").readline().rstrip()
        p = PersonalizeGenome(outDir, vcf, ref, hap1Ref, hap2Ref,rnaedit,editFile)
        p.personalize_genome()
        STAR_create_genome(outDir, ref, "REF",threads)
        STAR_create_genome(outDir, hap1Ref, "HAP1",threads)
        STAR_create_genome(outDir, hap2Ref, "HAP2",threads)



    elif setting == "mapping" :
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 3 :
        sys.stderr.write("Input arguments "+ command+" are not correct\n")
        sys.stderr.write(helpStr + '\n\n')
        sys.exit()
      elif len(command) == 3:
        if command[2].strip().lower() == "alleles":
          seqs = ".rPGASeqs.yaml"
          STAR_perform_mapping(outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
          STAR_perform_mapping(outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
          sam_to_sorted_bam(outDir+'/HAP1/STARalign/Aligned.out')
          sam_to_sorted_bam(outDir+'/HAP2/STARalign/Aligned.out')
        else:
          sys.stderr.write("Input arguments "+ command+" are not correct\n")
          sys.stderr.write(helpStr + '\n\n')
          sys.exit()
      else :
        seqs = ".rPGASeqs.yaml"
        STAR_perform_mapping(outDir, "REF", seqs,threads,mismatches,gzipped,multimapped)
        STAR_perform_mapping(outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
        STAR_perform_mapping(outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
        sam_to_sorted_bam(outDir+'/HAP1/STARalign/Aligned.out')
        sam_to_sorted_bam(outDir+'/HAP2/STARalign/Aligned.out')
        sam_to_sorted_bam(outDir+'/REF/STARalign/Aligned.out')
        

    elif setting == "discover" :
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2 :
        sys.stderr.write("Input arguments are not correct\n")
        sys.exit()
      else :
        discoverJunctions = True
        ref = open(".rPGAGenome.yaml").readline().rstrip()
        vcf = open(".rPGAGenotype.yaml").readline().rstrip()
        seqs = ".rPGASeqs.yaml"
        gtf = open(".rPGAJunctions.yaml").readline().rstrip()
        hap1Bam = outDir+'/HAP1/STARalign/Aligned.out.sorted.bam'
        hap2Bam = outDir+'/HAP2/STARalign/Aligned.out.sorted.bam'
        refBam = outDir+'/REF/STARalign/Aligned.out.sorted.bam'
        
        if multiprocessing:
          import multiprocessing
          CHROMS = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
          processes = [mp.Process(target=worker, args=(c,)) for c in CHROMS]
          for p in processes:
            p.start()
          for p in processes:
            p.join()

        else:
          p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting,rnaedit,editFile)
          p.haplotype_specific_junctions()

    elif setting == "alleles":
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2:
        sys.stderr.write("Input arguments are not correct\n")
        sys.exit()
      else:
        discoverJunctions = False        
        writeBam = True
        vcf = open(".rPGAGenotype.yaml").readline().rstrip()
        gtf = open(".rPGAJunctions.yaml").readline().rstrip()
        hap1Bam = outDir+'/HAP1/STARalign/Aligned.out.sorted.bam'
        hap2Bam = outDir+'/HAP2/STARalign/Aligned.out.sorted.bam'
        refBam = ""
        p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting,rnaedit,editFile)
        p.haplotype_specific_junctions()

    else :
      sys.stderr.write("rPGA run -- unknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")
    
