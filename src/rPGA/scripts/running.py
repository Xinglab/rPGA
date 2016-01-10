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
  p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, i, writeBam,discoverJunctions)
  p.haplotype_specific_junctions()

################################################################################
#             PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                 #
################################################################################

class PersonalizeGenome :
  def __init__(self, outDir, vcf, ref, hap1Ref, hap2Ref):
    self._outDir = outDir
    self._vcf = vcf
    self._ref = ref
    self._hap1Ref = hap1Ref
    self._hap2Ref = hap2Ref

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
                  v1[result['CHROM']][result['POS']] = [alleles[0],alleles[g1]]
                if (g2 != 0):
                  v2[result['CHROM']][result['POS']] = [alleles[0],alleles[g2]]
    vcf_in.close()
    return v1,v2

  def personalize_genome(self):
    hap1 = self.read_reference()
    vcf1,vcf2 = self.read_in_vcf()
    hap1_out = open(self._hap1Ref,'w')
    for chrom in vcf1:
      for pos in vcf1[chrom]:
        ref = vcf1[chrom][pos][0]
        alt = vcf1[chrom][pos][1]
        hap1['chr'+chrom][int(pos)-1] = alt
    for chrom in hap1:
      hap1_out.write('>'+chrom+'\n')
      hap1_out.write(''.join(hap1[chrom])+'\n')
    hap1.clear()
    hap1_out.close()

    hap2 = self.read_reference()
    hap2_out = open(self._hap2Ref,'w')
    for chrom in vcf2:
      for pos in vcf2[chrom]:
        ref = vcf2[chrom][pos][0]
        alt = vcf2[chrom][pos][1]
        hap2['chr'+chrom][int(pos)-1] = alt
    for chrom in hap2:
      hap2_out.write('>'+chrom+'\n')
      hap2_out.write(''.join(hap2[chrom])+'\n')
    hap2.clear()
    hap2_out.close()

    return

class DiscoverSpliceJunctions :
  def __init__(self, outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting):
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

  def read_in_vcf(self):
    vcf_in = open(self._vcf)
    v = defaultdict(list)
    vids = defaultdict(str)
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
          genopos += int(cigar[1]) # otherwise, continue searching throught the rest of the read
          readpos += int(cigar[1])
      elif cigar[0]==3: # intron                     
        genopos += int(cigar[1])
        if p <= (genopos): # snp falls into intron
          return 1
      else: #read has non N/M in cigar, throw it out                                              
        return 2
    return 3 # read is not allele specific

  def num_mismatches(self,r):
    # return number of mismatches
    nm = 0
    for tag in r.tags:
      if tag[0]=='NM': # NM is tag for mismatches
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
      return ','.join([v[i] for i in variant_positions])
    else:
      return []

  def haplotype_specific_junctions(self):
    bam1 = pysam.Samfile(self._hap1Bam)
    bam2 = pysam.Samfile(self._hap2Bam)
    hap1,hap2,conflicting = list(),list(),list()
    snpreads1,snpreads2 = defaultdict(lambda: defaultdict(list)), defaultdict(lambda: defaultdict(list))
    reads1,reads2 = defaultdict(lambda: defaultdict(list)),defaultdict(lambda: defaultdict(list))
    hetsnps,snpids = self.read_in_vcf()
    geneGroup,gtf,geneInfo = self.read_in_gtf()

    for p in sorted(hetsnps):
      for pileupcolumn in bam1.pileup('chr'+self._chromosome,int(p),int(p)+1):
        if int(pileupcolumn.pos)== p:
          for pileupread in pileupcolumn.pileups: 
            r = pileupread.alignment 
            snpreads1[r.pos][r.qname].append(p) 
            reads1[r.pos][r.qname] = r

      for pileupcolumn in bam2.pileup('chr'+self._chromosome,int(p),int(p)+1):
        if int(pileupcolumn.pos)== p:
          for pileupread in pileupcolumn.pileups:
            r = pileupread.alignment 
            snpreads2[r.pos][r.qname].append(p) 
            reads2[r.pos][r.qname] = r
            
    for pos in snpreads1:
      for qname in snpreads1[pos]: 
        if qname in snpreads2[pos]: 
          if self.is_unique(reads1[pos][qname]) and self.is_unique(reads2[pos][qname]):
            if snpreads1[pos][qname]==snpreads2[pos][qname]:
              check1,check2 = [],[] 
              for snp in snpreads1[pos][qname]:
                check1.append(self.check_cigar(reads1[pos][qname], hetsnps[snp][0], snp)) 
                check2.append(self.check_cigar(reads2[pos][qname], hetsnps[snp][1], snp)) 
              if ((1 not in check1+check2) and (2 not in check1+check2)):
                yes1 = sum([1 if int(i)==0 else 0 for i in check1])
                yes2 = sum([1 if int(i)==0 else 0 for i in check2])
                no1 = sum([1 if int(i)==3 else 0 for i in check1])
                no2 = sum([1 if int(i)==3 else 0 for i in check2]) 
                mismatch1 = self.num_mismatches(reads1[pos][qname]) 
                mismatch2 = self.num_mismatches(reads2[pos][qname])
                if ((yes1 > no1) and (mismatch1 < mismatch2)): 
                  hap1.append(qname)
                elif ((yes2 > no2) and (mismatch2 < mismatch1)): 
                  hap2.append(qname)
                else:
                  conflicting.append(qname)

    bamr = pysam.Samfile(self._refBam,"rb")
    inboth = [r for r in hap1 if r in hap2]
    hap1 = [r for r in hap1 if ((r not in inboth) and (r not in conflicting))]
    hap2 = [r for r in hap2 if ((r not in inboth) and (r not in conflicting))] 
    conflicting += inboth 
    counts = {}
    counts['hap1'] = len(set(hap1))/2
    counts['hap2'] = len(set(hap2))/2
    counts['conflicting'] = len(conflicting)/2
    
    if (self._discoverJunctions):
      junctions = defaultdict(lambda: defaultdict(set))
    if self._writeBam:
      out1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.bam','wb',template=bam1)
      out2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.bam','wb',template=bam2)
    if self._conflicting:
      conflict1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.conflicting.bam','wb',template=bam1)
      conflict2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.conflicting.bam','wb',template=bam2)
      
    for r in bam1.fetch('chr'+str(self._chromosome)):
      if self._writeBam:
        if r.qname in hap1:
          out1.write(r)
      if self._conflicting:
        if r.qname in conflicting:
          conflict1.write(r)
      if self._discoverJunctions:
        if ((r.qname not in conflicting) and (r.qname not in hap2)):
          juncs = self.get_junction_coordinates(r)
          for j in juncs:
            start,end = j
            junctions['1'][start,end].add(r.pos)
    
    for r in bam2.fetch('chr'+self._chromosome):
      if self._writeBam:
        if ( r.qname in hap2):
          out2.write(r)
      if self._conflicting:
        if r.qname in conflicting:
          conflict2.write(r)
      if self._discoverJunctions:
        if ((r.qname not in conflicting) and (r.qname not in hap1)):
          juncs = self.get_junction_coordinates(r)
          for j in juncs:
            start,end = j
            junctions['2'][start,end].add(r.pos)

    
    if self._discoverJunctions:
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
      spec['1'] = [j for j in nR['1'] if ((j not in nR['2']) and (j not in nR['R']) and nR['1'][j]>1 and len(self.get_splicesite_snp(start,end,snpids))>0)]
      spec['2'] = [j for j in nR['2'] if ((j not in nR['1']) and (j not in nR['R']) and nR['2'][j]>1 and len(self.get_splicesite_snp(start,end,snpids))>0)]
      spec['12'] = [j for j in nR['1'] if ((j in nR['2']) and (j not in nR['R']) and nR['1'][j]>1 and nR['2'][j]>1 and len(self.get_splicesite_snp(start,end,snpids))>0)]
      spec['R'] = [j for j in nR['R'] if ((j not in nR['1']) and (j not in nR['2'])and nR['R'][j]>1 and len(self.get_splicesite_snp(start,end,snpids))>0)]
      bed = defaultdict(list)
      for h in ['1','2','R']:
        for start,end in spec[h]:
          num_overlapping_reads = sum([nR[h][j] for j in nR[h] 
                                     if ( j[0]<end and j[1]>start and
                                          self.characterize_junctions(self._chromosome,j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
          if num_overlapping_reads > 0:
            freq = float(nR[h][start,end])/float(num_overlapping_reads)
          else:
            freq = 1
          n_or_r,strand = self.characterize_junctions(self._chromosome,start,end,geneGroup,gtf,geneInfo)
          bed[h].append('chr'+self._chromosome+' '+str(start) + ' ' + str(end) + ' J_'+str(counter)+'_'+n_or_r+'_'+snp+ ' ' + str(strand) + ' ' +str(freq))

      for start,end in spec['12']:
        num_overlapping_reads1 = sum([nR['1'][j] for j in nR['1'] if 
                                    ( j[0]<end and j[1]>start and 
                                      self.characterize_junctions(self._chromosome,j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
        num_overlapping_reads2 = sum([nR['2'][j] for j in nR['2'] if 
                                    ( j[0]<end and j[1]>start and 
                                      self.characterize_junctions(self._chromosome,j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
        if num_overlapping_reads1 == 0:
          freq1=1.0
        else:
          freq1 = float(nR['1'][start,end])/float(num_overlapping_reads1)
        if num_overlapping_reads2==0:
          freq2 = 1.0
        else:
          freq2 = float(nR['2'][start,end])/float(num_overlapping_reads2)
        n_or_r,strand = self.characterize_junctions(self._chromosome,start,end,geneGroup,gtf,geneInfo)
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
            "To run discover, where $ is your prompt:                    \n" +\
            "                                                            \n" +\
            "$ rPGA run discover -c chrom                                \n" 

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

  chromosome = args.c
  writeBam = args.b
  multiprocessing = args.p
  gzipped = args.g
  writeConflicting = args.conflict
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
        p = PersonalizeGenome(outDir, vcf, ref, hap1Ref, hap2Ref)
        p.personalize_genome()
        STAR_create_genome(outDir, ref, "REF",threads)
        STAR_create_genome(outDir, hap1Ref, "HAP1",threads)
        STAR_create_genome(outDir, hap2Ref, "HAP2",threads)



    elif setting == "mapping" :
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2 :
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
          p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting)
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
        hap1Bam = outDir+'/HAP1/STARalign/Aligned.out.sorted.bam'
        hap2Bam = outDir+'/HAP2/STARalign/Aligned.out.sorted.bam'
        p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting)
        p.allele_specific_assignment()

    else :
      sys.stderr.write("rPGA genomes -- unnknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")
    
