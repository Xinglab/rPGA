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

# imports
import sys, os, shutil,gzip
import subprocess
import re,logging,time,datetime,commands,argparse,random
from collections import defaultdict
import pysam, pybedtools
from pysam import Samfile
from pybedtools import BedTool
from  more_itertools import unique_everseen

## function and class definitions

## function and class definitions

def bam2sort(prefix):
  sys.stdout.write("sort and index: "+ prefix + '\n')
  bam_fn = prefix + '.bam'
  pysam.sort(bam_fn, prefix+'.sorted')
  pysam.index(prefix+'.sorted.bam')
  return

class DiscoverSpliceJunctions :
  def __init__(self, outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, rnaedit,editFile,gzipped):
    self._outDir = outDir
    if os.path.isdir(vcf):
      self.vcf = vcf
      self.vcfdir = True
    else:
      self.vcf = vcf
      self.vcfdir = False
    self.gtf = gtf
    self.hap1Bam = hap1Bam
    self.hap2Bam = hap2Bam
    self.refBam = refBam
    self.VCFheader = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    self.rnaedit = rnaedit
    self.editFile = editFile
    self.gzipped = gzipped

  def get_tags(self,r): ## get read tags
    return {key: value for key, value in r.tags}

  def read_in_rna_editing(self):
    logging.info('reading in rna editing events')
    e = defaultdict(lambda:defaultdict(list))
    if not self.rnaedit:
      return e
    fin = open(self.editFile)
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
        e[result['chromosome'][3:]][(int(result['position'])-1)/100].append(int(result['position'])-1)
    return e

  def read_in_vcf(self):
    v = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    vids = defaultdict(lambda:defaultdict(lambda:defaultdict(str)))
    if self.vcfdir:
      fList = [os.path.join(self.vcf,f) for f in os.listdir(self.vcf)]
    else:
      fList = [self.vcf]
    for f in fList:
      if self.gzipped:
        vcf_in = gzip.open(f)
      else:
        vcf_in = open(f)
      sys.stdout.write("Reading in VCF: " + str(f) + "\n")
      for line in vcf_in:
        if line.startswith('#'): # skip over header lines
          continue
        result = {} # store line in dictionary
        fields = line.rstrip().split()
        for i,col in enumerate(self.VCFheader):
          result[col] = fields[i]
        infos = [x for x in result['INFO'].split(';') if x.strip()]
        for i in infos:
          if '=' in i:
            key,value = i.split('=') 
            result[key] = value
        alt = result['ALT'].split(',')
        alleles = [result['REF']] + alt
        if all([re.match(r'[ACGT]',i) for i in alleles]) and all([len(i)==1 for i in alleles]):
          geno = result['SAMPLE'].split(':')[0]
          if '|' in geno: ## if genotype is phased
            g1 = int(geno.split('|')[0]) ## hap1 genotype
            g2 = int(geno.split('|')[1]) ## hap2 genotype
            if (g1!=g2): # heterozygous snp
              v[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = [alleles[g1],alleles[g2],alleles[0]] ## subtract one from position (1 based) to match bam file (0 based)
              vids[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
            else:
              vids[result['CHROM']][(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
      vcf_in.close()
    return v,vids


  def read_in_gtf(self):
    gtf_in = open(self.gtf)
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

  def get_splicesite_snp(self,s,e,v):
    # returns splice site snp in junction start-end
    variant_positions = [p for p in (range(s-1,s+1)+range(e-2,e)) if p in v[p/100]]
    if len(variant_positions)>0:
      return [v[i/100][i] for i in variant_positions]
    else:
      return []

  def haplotype_specific_read(self,r,v,h,e):
    # return 0=hap spec, 1=not hap spec, 2=conflicting, 5 = doesn't cover het snp 
    reference_positions = r.get_reference_positions()
    edit = [p for p in reference_positions if p in e[p/100]]
    snps = [p for p in reference_positions if ((p in v[p/100]) and (p not in edit))]
    if len(snps)==0: # read doesn't cover het snp
      return 5, edit, [], []
    check = [0 if r.seq[reference_positions.index(p)]==v[p/100][p][h] else 3 for p in snps]
    allele = ['R' if v[p/100][p][h]==v[p/100][p][2] else 'A' for p in snps]
    if check.count(0)>check.count(3):  # majority votes for read being hap spec
      return 0,edit,snps,allele
    elif check.count(0)==check.count(3): # conflicting
      return 2,edit,snps,allele
    else: #majority votes for other haplotype
      return 1,edit,snps,allele




  def discover_junctions_bychrom(self,chrom,hetsnps,snpids,editpositions):
    bam1 = pysam.Samfile(self.hap1Bam)
    bam2 = pysam.Samfile(self.hap2Bam)
    snpreads1,snpreads2 = defaultdict(list),defaultdict(list)
    snpreads = defaultdict(lambda:defaultdict(list))
    reads1,reads2 = defaultdict(list),defaultdict(list)

    spec1,spec2  = list(),list()
    snps1,snps2 = defaultdict(list),defaultdict(list)
    refalt1,refalt2 = defaultdict(list),defaultdict(list)
    edit1,edit2= defaultdict(list),defaultdict(list)
    hap1only,hap2only = list(),list()
    sys.stdout.write( "Reading in hap1 bam file: "+ self.hap1Bam + "\n")
    
    for r in bam1.fetch('chr'+str(chrom)):
      tags = self.get_tags(r)
      if int(tags['NH'])>1: ## read is multimapped, deal with separately
        continue
      spec,editpos,snppos,refalt = self.haplotype_specific_read(r,hetsnps,0, editpositions)
      snpreads1[spec].append(r.qname)
      reads1[r.qname].append(r)
      edit1[r.qname] += editpos
      snps1[r.qname] += snppos
      refalt1[r.qname] += refalt

    sys.stdout.write( "Reading in hap2 bam file: "+ self.hap2Bam+"\n")
    for r in bam2.fetch('chr'+str(chrom)):
      tags = self.get_tags(r)
      if int(tags['NH'])>1: ## read is multimapped
        continue
      spec,editpos,snppos,refalt = self.haplotype_specific_read(r,hetsnps,1, editpositions)
      snpreads2[spec].append(r.qname)
      reads2[r.qname].append(r)
      edit2[r.qname]+=editpos
      snps2[r.qname]+=snppos
      refalt2[r.qname] += refalt

    sys.stdout.write("Assign haplotype specific reads\n")
    conflicting = list(set([ r for r in snpreads1[0] if r in snpreads2[0]] + snpreads1[2] + snpreads2[2]))
    multimapped = list(set(snpreads1[3] + snpreads2[3]))
    rnaeditreads = list(set(snpreads1[4] + snpreads2[4]))
    snpreads1[0] = list(set([ r for r in snpreads1[0] if ((r not in conflicting) and (r not in snpreads2[3]))]))
    snpreads2[0] = list(set([ r for r in snpreads2[0] if ((r not in conflicting) and (r not in snpreads1[3]))]))

    for qname in snpreads1[0]:
      if qname in reads2:
        if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads1[qname]))]): # reads have same starting position in hap1 and hap2
          if self.num_mismatches(reads1[qname][0]) < self.num_mismatches(reads2[qname][0]):
            spec1.append(qname)
      else:
        hap1only.append(qname)
      
    for qname in snpreads2[0]:
      if qname in reads1:
        if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads2[qname]))]): # reads have same starting position in hap1 and hap2 
          if self.num_mismatches(reads2[qname][0]) < self.num_mismatches(reads1[qname][0]):
            spec2.append(qname)
      else:
        hap2only.append(qname)

    conflictCount = len(set(conflicting))
    hap1Count = len(set(spec1))
    hap2Count = len(set(spec2))

    geneGroup,gtf,geneInfo = self.read_in_gtf()
    bamr = pysam.Samfile(self.refBam,"rb")
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

    for r in bamr.fetch('chr'+str(chrom)):
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
                                          self.characterize_junction(str(chrom),j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
        if num_overlapping_reads > 0:
          freq = float(nR[h][start,end])/float(num_overlapping_reads)
        else:
          freq = 1
        n_or_r,strand = self.characterize_junction(str(chrom),start,end,geneGroup,gtf,geneInfo)
        snp = ','.join(self.get_splicesite_snp(start,end,snpids))
        bed[h].append('chr'+str(chrom)+' '+str(start) + ' ' + str(end) + ' J_'+str(counter)+'_'+n_or_r+'_'+snp+ ' ' + str(strand) + ' ' +str(freq))
          
    counter = 0
    for start,end in spec['12']:
      counter += 1
      num_overlapping_reads1 = sum([nR['1'][j] for j in nR['1'] if 
                                    ( j[0]<end and j[1]>start and 
                                      self.characterize_junction(str(chrom),j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
      num_overlapping_reads2 = sum([nR['2'][j] for j in nR['2'] if 
                                    ( j[0]<end and j[1]>start and 
                                      self.characterize_junction(str(chrom),j[0],j[1],geneGroup,gtf,geneInfo)[0]=="R")])
      if num_overlapping_reads1 == 0:
        freq1=1.0
      else:
        freq1 = float(nR['1'][start,end])/float(num_overlapping_reads1)
      if num_overlapping_reads2==0:
        freq2 = 1.0
      else:
        freq2 = float(nR['2'][start,end])/float(num_overlapping_reads2)
      n_or_r,strand = self.characterize_junction(str(chrom),start,end,geneGroup,gtf,geneInfo)
      snp =','.join(self.get_splicesite_snp(start,end,snpids))
      bed['12'].append('chr'+str(chrom)+' '+str(start) + ' ' + str(end) + ' J_'+str(counter)+'_'+n_or_r+'_'+snp+ ' ' + str(strand)+' '+str((freq1+freq2)/2))

    return '\n'.join(bed['1']),'\n'.join(bed['2']),'\n'.join(bed['12']), '\n'.join(bed['R'])



  def discover_junctions(self):
    sys.stdout.write("Reading in VCF and editing files\n")
    VCF,VCFids = self.read_in_vcf()
    RNAedit = self.read_in_rna_editing()
    bed1, bed2, bed12, bedR = defaultdict(str),defaultdict(str),defaultdict(str),defaultdict(str)
    for c in VCF:
      sys.stdout.write("Discovering junctions in chromosome "+str(c) + '\n')
      bed1[c], bed2[c],bed12[c],bedR[c] = self.discover_junctions_bychrom(c,VCF[c],VCFids[c],RNAedit[c])
    bedstring1, bedstring2, bedstring12, bedstringR = "", "", "", ""
    for c in VCF:
      bedstring1 += (bed1[c] + '\n')
      bedstring2 += (bed2[c] + '\n')
      bedstring12 += (bed12[c] + '\n')
      bedstringR += (bedR[c] + '\n')
    mybed1 = BedTool(bedstring1, from_string=True).sort().saveas(self._outDir + '/hap1.specific.bed')
    mybed2 = BedTool(bedstring2, from_string=True).sort().saveas(self._outDir + '/hap2.specific.bed')
    mybed12 = BedTool(bedstring12, from_string=True).sort().saveas(self._outDir + '/hap1hap2.specific.bed')
    mybedR = BedTool(bedstringR, from_string=True).sort().saveas(self._outDir + '/ref.specific.bed')
    return


def main(args):

  outDir = args.o
  vcf = args.v
  if args.b1:
    hap1Bam = args.b1
  else:
    hap1Bam = os.path.join(outDir,"HAP1","STARalign","Aligned.out.sorted.bam")
  if args.b2:
    hap2Bam = args.b2
  else:
    hap2Bam = os.path.join(outDir,"HAP2","STARalign","Aligned.out.sorted.bam")
  if args.br:
    refBam = args.br
  else:
    refBam = os.path.join(outDir,"REF","STARalign","Aligned.out.sorted.bam")

  rnaedit = args.rnaedit
  editFile = args.e
  gzipped = args.gz
  gtf = args.g

  if not args.o:
    sys.stderr.write("rPGA: ERROR -o outdir option is required!\n\n")
    sys.exit()
  if not args.v:
    sys.stderr.write("rPGA: ERROR -v vcf option is required!\n\n")
    sys.exit()
  if not args.g:
    sys.stderr.write("rPGA: ERROR -v gtf option is required!\n\n")
    sys.exit()

  p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, rnaedit,editFile,gzipped)
  p.discover_junctions()
