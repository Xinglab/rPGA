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
import sys, os, shutil,gzip
import progressbar
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage,  ProgressBar, ReverseBar, RotatingMarker,SimpleProgress, Timer
import subprocess
import re,logging,time,datetime,commands,argparse,random
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
  opts += (seqs + ' ')
#  for line in open(seqs) :
#    line = line.rstrip()
#    opts += (str(line) + " ")
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
  def __init__(self, outDir, vcf, ref, hap1Ref, hap2Ref,rnaedit,editFile,gzipped):
    self._outDir = outDir
    self._vcf = vcf
    self._ref = ref
    self._hap1Ref = hap1Ref
    self._hap2Ref = hap2Ref
    self._rnaedit = rnaedit
    self._editFile = editFile
    self._report = os.path.join(outDir,'report.personalize.txt')
    self._widgets = [Percentage(),' Processed: ', Counter(), ' lines (', Timer(), ')']
    self._gzipped = gzipped
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
    print " read in edit file"
    edit_in = open(self._editFile)
    sys.stdout.write('# reading rna-edit file: '+self._editFile+'\n')
#    num_lines = sum(1 for line in open(self._editFile))
#    print num_lines
#    pbar = ProgressBar(widgets=self._widgets,max_value=num_lines)
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
    print "# storing vcf files"
    sys.stdout.write('# storing vcf files \n')
    CHROMS = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
#    CHROMS = ['21']
    VCF_HEADER = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    v1,v2 = defaultdict(lambda: defaultdict(list)),defaultdict(lambda: defaultdict(list))
    for c in CHROMS:
      if self._gzipped:
        vcf_fn = self._vcf + '/' + c + '.vcf.gz'
        vcf_in = gzip.open(vcf_fn)
      else:
        vcf_fn = self._vcf + '/' + c +'.vcf'
        vcf_in = open(vcf_fn)
      sys.stdout.write('# reading in '+vcf_fn+'\n') 
#      num_lines = sum(1 for line in open(vcf_fn))
#      pbar = ProgressBar(widgets=self._widgets,max_value=num_lines)
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
    print "personalizing genome"
    if self._rnaedit:
      hap1 = self.read_reference()
      vcf1,vcf2 = self.read_in_vcf()
      editSites = self.read_in_edit()
      report_out = open(self._report,'w')
      hap1_out = open(self._hap1Ref,'w')
      editCounter = 0
      snpCounter = 0
      edit_snp_sites = []
      for chrom in vcf1:
        for pos in vcf1[chrom]:
          if pos in editSites[chrom]:
            editCounter += 1
            hap1['chr'+chrom][int(pos)-1] = 'N'
            edit_snp_sites.append('chr'+chrom+' '+str(pos))
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
      report_out.write('# number of hap1 SNPs not overlapping RNA editing sites: ' + str(snpCounter) + '\n')
      report_out.write('# total number of changed bases in hap1: '+str(editCounter + snpCounter) + '\n')
      report_out.write('# overlapping editing sites and overlapping SNPs in hap1:\n')
      report_out.write('\n'.join(edit_snp_sites))

      editCounter = 0
      snpCounter = 0
      hap2 = self.read_reference()
      hap2_out = open(self._hap2Ref,'w')
      edit_snp_sites = []
      for chrom in vcf2:
        for pos in vcf2[chrom]:
          if pos in editSites[chrom]:
            editCounter+=1
            hap2['chr'+chrom][int(pos)-1] = 'N'
            edit_snp_sites.append('chr'+chrom+' '+str(pos))
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
      report_out.write('# number of hap2 SNPs not overlapping RNA editing sites: ' + str(snpCounter) + '\n')
      report_out.write('# total number of changed bases in hap2: '+str(editCounter + snpCounter) + '\n')
      report_out.write('# overlapping editing sites and overlapping SNPs in hap2:\n')
      report_out.write('\n'.join(edit_snp_sites))
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
  def __init__(self, outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting,rnaedit,editFile,gzipped,printall,consensus):
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
    self._gzipped = gzipped
    self._printall = printall
    self._consensus = consensus
    self._widgets = [Percentage(),' Processed: ', Counter(), ' lines (', Timer(), ')']
    
  def read_in_rna_editing(self):
    logging.info('reading in rna editing events')
    e = defaultdict(lambda:defaultdict(list))
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
        e[result['chromosome'][3:]][(int(result['position'])-1)/100].append(int(result['position'])-1)
    return e

  def read_in_vcf(self):
    if self._gzipped:
      vcf_in = gzip.open(self._vcf+'.gz')
#      num_lines = sum(1 for line in gzip.open(self._vcf + '.gz'))
    else:
      vcf_in = open(self._vcf)
#      num_lines = sum(1 for line in open(self._vcf)) 
#    pbar = ProgressBar(widgets=self._widgets,max_value=num_lines) 
    v = defaultdict(lambda:defaultdict(list))
    vids = defaultdict(lambda:defaultdict(str))

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
                    v[(int(result['POS'])-1)/100][int(result['POS'])-1] = [alleles[g1],alleles[g2],alleles[0]] ## subtract one from position (1 based) to match bam file (0 based)
                    vids[(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
                  else:
                    vids[(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']


    else:
      editPos = defaultdict(lambda:defaultdit(list))
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
                  v[(int(result['POS'])-1)/100][int(result['POS'])-1] = [alleles[g1],alleles[g2],alleles[0]] ## subtract one from position (1 based) to match bam file (0 based)                                                                 
                  vids[(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']
                else:
                  vids[(int(result['POS'])-1)/100][int(result['POS'])-1] = result['ID']

    vcf_in.close()
    return v,vids,editPos
    
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

  def get_splicesite_snp(self,s,e,v):
    # returns splice site snp in junction start-end
    variant_positions = [p for p in (range(s-1,s+1)+range(e-2,e)) if p in v[p/100]]
    if len(variant_positions)>0:
      return [v[i/100][i] for i in variant_positions]
    else:
      return []

  def haplotype_specific_read(self,r,v,h,e):
    # return 0 = hap spec, 1 = not hap spec, 2 = conflicting, 3 = multimapped, 4 = rnaedit,  5 = doesn't cover het snp
    if not self.is_unique(r):
      return 3,[],[]
    reference_positions = r.get_reference_positions()
    snps = [p for p in reference_positions if p in v[p/100]]
    edit = [p for p in reference_positions if p in e[p/100]]
    if len(edit)>0: # read covers rna-editing position
      return 4,edit,[]
    elif len(snps)==0: # read doesn't cover het snp
      return 5,[],[]
    else:
      check = [0 if r.seq[reference_positions.index(p)]==v[p/100][p][h] else 3 for p in snps]
#      check = [self.check_cigar(r, v[p/100][p][h], p) for p in snps]
      allele = ['R' if v[p/100][p][h]==v[p/100][p][2] else 'A' for p in snps]
      if check.count(0)>check.count(3):
        return 0,snps,allele
      elif check.count(0)==check.count(3):
        return 2,snps,allele
      else:
        return 1,snps,allele

      

  def haplotype_specific_junctions(self):
    print "assigning haplotype specific reads"
    bam1 = pysam.Samfile(self._hap1Bam)
    bam2 = pysam.Samfile(self._hap2Bam)
    snpreads1,snpreads2 = defaultdict(list),defaultdict(list)
    snpreads = defaultdict(lambda:defaultdict(list))
    reads1,reads2 = defaultdict(list),defaultdict(list)

    print "## reading in VCF file"
    hetsnps,snpids,editpositions = self.read_in_vcf()
    spec1,spec2  = list(),list()
    snps1,snps2 = defaultdict(list),defaultdict(list)
    refalt1,refalt2 = defaultdict(list),defaultdict(list)
    hap1only,hap2only = list(),list()
    consensus_reads = defaultdict(lambda:defaultdict(list))

    print "## reading in hap1 bam file: ", self._hap1Bam
    for r in bam1.fetch('chr'+str(self._chromosome)):
      spec,snppos,refalt = self.haplotype_specific_read(r,hetsnps,0, editpositions[self._chromosome])
      snpreads1[spec].append(r.qname)
      reads1[r.qname].append(r)
      snps1[r.qname] += snppos
      refalt1[r.qname] += refalt
      consensus_reads[r.qname][1].append(spec)

    print "## reading in hap2 bam file: ", self._hap2Bam
    for r in bam2.fetch('chr'+str(self._chromosome)):
      spec,snppos,refalt = self.haplotype_specific_read(r,hetsnps,1, editpositions[self._chromosome])
      snpreads2[spec].append(r.qname)
      reads2[r.qname].append(r)
      snps2[r.qname]+=snppos
      refalt2[r.qname] += refalt
      consensus_reads[r.qname][2].append(spec)

    conflicting = list(set([ r for r in snpreads1[0] if r in snpreads2[0]] + snpreads1[2] + snpreads2[2]))
    multimapped = list(set(snpreads1[3] + snpreads2[3]))
    rnaeditreads = list(set(snpreads1[4] + snpreads2[4]))
    snpreads1[0] = list(set([ r for r in snpreads1[0] if ((r not in conflicting) and (r not in snpreads2[3]))]))
    snpreads2[0] = list(set([ r for r in snpreads2[0] if ((r not in conflicting) and (r not in snpreads1[3]))]))
    

    print "## assigning specific reads"
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

    print "## printing report file"
    report_out = open(self._outDir + '/report.'+self._chromosome+'.txt','w')
    report_out.write('########## haplotype assignment report for chromosome ' + self._chromosome + '##########\n')
    report_out.write('# hap1' + '\t' + str(hap1Count) + '\n')
    report_out.write('# hap2' +'\t' + str(hap2Count) + '\n')
    report_out.write('# conflicting' +'\t' + str(conflictCount) + '\n')
    
    if self._rnaedit:
      print "printing rna editing file"
      edit1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.rnaedit.bam','wb',template=bam1)
      edit2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.rnaedit.bam','wb',template=bam2)
      report_out.write('# rna-edit,hap1\t' + str(len(set(snpreads1[4]))) + '\n')
      report_out.write('# rna-edit,hap2\t' + str(len(set(snpreads2[4]))) + '\n')
      for qname in snpreads1[4]:
        for r in reads1[qname]:
          r.tags += [('HT',1),('EP',';'.join([str(s) for s in snps1[qname]]))]
          edit1.write(r)
      for qname in snpreads2[4]:
        for r in reads2[qname]:
          r.tags += [('HT',2),('EP',';'.join([str(s) for s in snps2[qname]]))]
          edit2.write(r)
    report_out.write('###################################################################################\n')
    if self._writeBam:
      print "## printing bam file"
      out1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.bam','wb',template=bam1)
      out2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.bam','wb',template=bam2)
      for qname in spec1:
        for r in reads1[qname]:
          r.tags += [('HT',1),('SP',';'.join([str(p) for p in snps1[qname]])),('GT',';'.join([x for x in refalt1[qname]]))]
          out1.write(r)

      for qname in spec2:
        for r in reads2[qname]:
          r.tags += [('HT',2),('SP',';'.join([str(p) for p in snps2[qname]])),('GT',';'.join([x for x in refalt2[qname]]))]
          out2.write(r)
      
      if self._printall: # print all reads that dont cover het snps
        for qname in snpreads1[5]:
          for r in reads1[qname]:
            r.tags += [('HT',1)]
            out1.write(r)

        for qname in snpreads2[5]:
          for r in reads2[qname]:
            r.tags += [('HT',2)]
            out2.write(r)

    if self._consensus: # print one consensus bam file
      print "## printing consensus bam file"
      multi = len(multimapped)
      rna_edit = len(rnaeditreads)
      hap1_only = 0
      hap2_only = 0
      dif_loc = 0
      hap1_spec = 0
      hap2_spec = 0
      conflict = 0
      same_hap1 = 0
      same_hap2 = 0
      consensus = pysam.Samfile(self._outDir + "/consensus."+self._chromosome+'.bam','wb',template=bam1)
      #print hap1 specific reads
      print "print hap1 specific reads to consensus file"
      if self._writeBam:
        for qname in spec1:
          hap1_spec += 1
          for r in reads1[qname]:
            consensus.write(r)
      else:
        for qname in spec1:
          hap1_spec += 1
          for r in reads1[qname]:
            r.tags += [('HT',1),('SP',';'.join([str(p) for p in snps1[qname]])),('GT',';'.join([x for x in refalt1[qname]]))]
            consensus.write(r)
      for qname in hap1only:
        hap1_only += 1
        for r in reads1[qname]:
          if len(snps1[qname])==0:
            r.tags += [('HT',1)]
          else:
            r.tags += [('HT',1),('SP',';'.join([str(p) for p in snps1[qname]])),('GT',';'.join([x for x in refalt1[qname]]))]
          consensus.write(r)
      #print hap2 specific reads
      print "print hap2 specific reads to consensus file"
      if self._writeBam:
        for qname in spec2:
          hap2_spec+=1
          for r in reads2[qname]:
           consensus.write(r)
      else:
        for qname in spec2:
          hap2_spec += 1
          for r in reads2[qname]:
            r.tags += [('HT',2),('SP',';'.join([str(p) for p in snps2[qname]])),('GT',';'.join([x for x in refalt2[qname]]))]
            consensus.write(r)
      for qname in hap2only:
        hap2_only += 1
        if len(snps2[qname])==0:
          r.tags += [('HT',2)]
        else:
          r.tags += [('HT',2),('SP',';'.join([str(p) for p in snps2[qname]])),('GT',';'.join([x for x in refalt2[qname]]))]
        consensus.write(r)
      #print conflicting reads (choose which haplotype alignment to print randomly)
      print "print conflicting reads to consensus file"
      for qname in conflicting:
        conflict += 1
        x = random.random()
        if x<0.5:
          for r in reads1[qname]:
            r.tags += [('HT',1),('SP',';'.join([str(p) for p in snps1[qname]])),('GT',';'.join([x for x in refalt1[qname]]))]
            consensus.write(r)
        else:
          for r in reads2[qname]:
            r.tags += [('HT',2),('SP',';'.join([str(p) for p in snps2[qname]])),('GT',';'.join([x for x in refalt2[qname]]))]
            consensus.write(r)

      #print reads that dont cover heterozygous snp(s)
      print "reads that dont cover het snps to consensus file"
      for qname in consensus_reads:
        if qname not in multimapped and qname not in rnaeditreads:
          if 1 in consensus_reads[qname]:
            if 2 in consensus_reads[qname]:
              if consensus_reads[qname][1]==[5,5] and consensus_reads[qname][2]==[5,5]: #both reads are uniquely mapped and dont cover het snp
                if all([True if reads1[qname][i].pos==reads2[qname][i].pos else False for i in range(len(reads2[qname]))]): #reads have same start
                  x = random.random()
                  if x < 0.5:
                    same_hap1 += 1
                    for r in reads1[qname]:
                      r.tags += [('HT',1)]
                      consensus.write(r)
                  else:
                    same_hap2 += 1
                    for r in reads2[qname]:
                      r.tags += [('HT',2)]
                      consensus.write(r)
                else:
                  dif_loc += 1
            else: #read maps to hap1 only
              hap1_only += 1
              for r in reads1[qname]:
                r.tags += [('HT',1)]
                consensus.write(r)
          else:
            hap2_only += 1
            for r in reads2[qname]:
              r.tags += [('HT',1)]
              consensus.write(r)
      report_out.write('####### consensus BAM file report #######\n')
      report_out.write('# total read pairs in consensus.'+self._chromosome+'.bam: '+str(hap1_only+hap2_only+hap1_spec+hap2_spec+conflict+same_hap1+same_hap2) + '\n')
      report_out.write('# - only mapped to hap1: '+str(hap1_only)+'\n')
      report_out.write('# - only mapped to hap2: '+str(hap2_only)+'\n')
      report_out.write('# - hap1 specific: ' + str(hap1_spec) + '\n')
      report_out.write('# - hap2 specific: ' + str(hap2_spec) + '\n')
      report_out.write('# - conflicting: ' + str(conflict) + '\n')
      report_out.write('# - same alignment, choose hap1: ' + str(same_hap1) + '\n')
      report_out.write('# - same alignment, choose hap2: ' + str(same_hap2) + '\n')
      report_out.write('# total read pairs not mapped to consensus.'+self._chromosome+'.bam: ' + str(multi+dif_loc+rna_edit) + '\n')
      report_out.write('# - multimapped: ' + str(multi) + '\n')
      report_out.write('# - mapped to different locations in each haplotype: ' + str(dif_loc) + '\n')
      report_out.write('# - cover rna-editing site: ' + str(rna_edit) + '\n')
      report_out.write('##########################################\n')

    if self._conflicting:
      print "## printing conflicting file"
      conflict1 = pysam.Samfile(self._outDir + "/hap1."+self._chromosome+'.conflicting.bam','wb',template=bam1)
      conflict2 = pysam.Samfile(self._outDir + "/hap2."+self._chromosome+'.conflicting.bam','wb',template=bam2)
      for qname in conflicting:
        for r in reads1[qname]:
          conflict1.write(r)
        for r in reads2[qname]:
          conflict2.write(r)
    
    if self._discoverJunctions:
      geneGroup,gtf,geneInfo = self.read_in_gtf()
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

  printall = args.printall
  writeBam = args.writeBam
#  multiprocessing = args.p
  gzipped = args.gz
  writeConflicting = args.conflict
  rnaedit = args.rnaedit
  consensus = args.consensus
  outDir = args.o

  if not os.path.exists(outDir):
    os.makedirs(outDir)

  ref = args.r
  vcf = args.v
  gtf = args.g

  if len(command) == 1 or (len(command)==2 and isHelpString(command[1].strip().lower())):
    sys.stderr.write(helpStr + "\n\n")
  else :
    setting = command[1].strip().lower()
#    outDir = open(".rPGAProject.yaml").readline().rstrip()
    hap1Ref = os.path.join(outDir, "hap1.fa")
    hap2Ref = os.path.join(outDir, "hap2.fa")

    if setting == "personalize" :
      print "personalizing genome"
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2 :
        sys.stderr.write("ERROR: Input arguments are not correct\n")
        sys.stderr.write(helpStr + '\n\n')
        sys.exit()
      else :
        if not args.r:
          sys.stderr.write("ERROR: rPGA run personalize command requires -r parameter \nExample: rPGA run mapping -v vcf_directory -r reference.fa  -o rPGA \n")
          sys.exit()
        if not args.v:
          sys.stderr.write("ERROR: rPGA run personalize command requires -v parameter \nExample: rPGA run mapping -v vcf_directory -r reference.fa  -o rPGA \n")
          sys.exit()
#        ref = open(".rPGAGenome.yaml").readline().rstrip()
#        vcf = open(".rPGAGenotype.yaml").readline().rstrip()
        p = PersonalizeGenome(outDir, vcf, ref, hap1Ref, hap2Ref,rnaedit,editFile,gzipped)
        p.personalize_genome()


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
#          seqs = ".rPGASeqs.yaml"
          seqs = ' '.join((args.s).split(','))
          if len((args.s).split(','))==0 or len((args.s).split(','))>2:
            sys.stderr.write("ERROR: Sequence parameter -s input is  not correct\n Example: rPGA run mappng alleles -s reads_1.fq,reads_2.fq -o rPGA\n")
            sys.exit()
#          ref = open(".rPGAGenome.yaml").readline().rstrip()
#          vcf = open(".rPGAGenotype.yaml").readline().rstrip()
          if not os.path.exists(os.path.join(outDir, "HAP1/STARindex")):
            os.makedirs(os.path.join(outDir, "HAP1/STARindex"))
          if not os.path.exists(os.path.join(outDir, "HAP2/STARindex")):
            os.makedirs(os.path.join(outDir, "HAP2/STARindex"))
          if not os.path.exists(os.path.join(outDir, "HAP1/STARalign")):
            os.makedirs(os.path.join(outDir, "HAP1/STARalign"))
          if not os.path.exists(os.path.join(outDir, "HAP2/STARalign")):
            os.makedirs(os.path.join(outDir, "HAP2/STARalign"))
          STAR_create_genome(outDir, hap1Ref, "HAP1",threads)
          STAR_create_genome(outDir, hap2Ref, "HAP2",threads)
          STAR_perform_mapping(outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
          STAR_perform_mapping(outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
          sam_to_sorted_bam(outDir+'/HAP1/STARalign/Aligned.out')
          sam_to_sorted_bam(outDir+'/HAP2/STARalign/Aligned.out')
          os.remove(os.path.join(outDir,'HAP1/STARalign/Aligned.out.sam'))
          os.remove(os.path.join(outDir,'HAP2/STARalign/Aligned.out.sam'))
          os.remove(os.path.join(outDir,'HAP1/STARalign/Aligned.out.bam'))
          os.remove(os.path.join(outDir,'HAP2/STARalign/Aligned.out.bam'))
        else:
          sys.stderr.write("Input arguments "+ command+" are not correct\n")
          sys.stderr.write(helpStr + '\n\n')
          sys.exit()
      else :
#        seqs = ".rPGASeqs.yaml"
        if not args.r:
          sys.stderr.write("ERROR: rPGA run mapping command requires -r parameter \nExample: rPGA run mapping -r reference.fa -s reads_1.fq,reads_.fq -o rPGA \n")
          sys.exit()
        seqs = ' '.join((args.s).split(','))
        if len((args.s).split(','))==0 or len((args.s).split(','))>2:
          sys.stderr.write("ERROR: Sequence parameter -s input is  not correct\n Example: rPGA run mapping - r reference.fa -s reads_1.fq,reads_2.fq -o rPGA\n")
          sys.exit()
        print "checking STAR index directories"
        if not os.path.exists(os.path.join(outDir, "HAP1/STARindex")):
          os.makedirs(os.path.join(outDir, "HAP1/STARindex"))
        if not os.path.exists(os.path.join(outDir, "HAP2/STARindex")):
          os.makedirs(os.path.join(outDir, "HAP2/STARindex"))
        if not os.path.exists(os.path.join(outDir, "REF/STARindex")):
          os.makedirs(os.path.join(outDir, "REF/STARindex"))
        if not os.path.exists(os.path.join(outDir, "HAP1/STARalign")):
          os.makedirs(os.path.join(outDir, "HAP1/STARalign"))
        if not os.path.exists(os.path.join(outDir, "HAP2/STARalign")):
          os.makedirs(os.path.join(outDir, "HAP2/STARalign"))
        if not os.path.exists(os.path.join(outDir, "REF/STARalign")):
          os.makedirs(os.path.join(outDir, "REF/STARalign"))
        print "creating STAR genome indicies"
        STAR_create_genome(outDir, ref, "REF",threads)
        STAR_create_genome(outDir, hap1Ref, "HAP1",threads)
        STAR_create_genome(outDir, hap2Ref, "HAP2",threads)        
        print "perform STAR mapping"
        STAR_perform_mapping(outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
        STAR_perform_mapping(outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
        STAR_perform_mapping(outDir, "REF", seqs,threads,mismatches,gzipped,multimapped)
        sam_to_sorted_bam(outDir+'/HAP1/STARalign/Aligned.out')
        sam_to_sorted_bam(outDir+'/HAP2/STARalign/Aligned.out')
        sam_to_sorted_bam(outDir+'/REF/STARalign/Aligned.out')
        os.remove(os.path.join(outDir,'HAP1/STARalign/Aligned.out.sam'))
        os.remove(os.path.join(outDir,'HAP2/STARalign/Aligned.out.sam'))
        os.remove(os.path.join(outDir,'REF/STARalign/Aligned.out.sam'))
        os.remove(os.path.join(outDir,'HAP1/STARalign/Aligned.out.bam'))
        os.remove(os.path.join(outDir,'HAP2/STARalign/Aligned.out.bam'))
        os.remove(os.path.join(outDir,'REF/STARalign/Aligned.out.bam'))

    elif setting == "discover" :
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2 :
        sys.stderr.write("ERROR: Input arguments are not correct\n")
        sys.exit()
      else :
        discoverJunctions = True
#        ref = open(".rPGAGenome.yaml").readline().rstrip()
#        vcf = open(".rPGAGenotype.yaml").readline().rstrip()
#        seqs = ".rPGASeqs.yaml"
#        gtf = open(".rPGAJunctions.yaml").readline().rstrip()
        if not args.v:
          sys.stderr.write("ERROR: rPGA run discover command requires -v parameter \nExample: rPGA run discover -c 1 -v vcf_directory -g annotation.gtf  -o rPGA \n")
          sys.exit()
        if not args.g:
          sys.stderr.write("ERROR: rPGA run discover command requires -g parameter \nExample: rPGA run discover -c 1 -v vcf_directory -g annotation.gtf  -o rPGA \n")
          sys.exit()
        if args.b1:
          hap1Bam = args.b1
        else:
          hap1Bam = outDir+'/HAP1/STARalign/Aligned.out.sorted.bam'
        if args.b2:
          hap2Bam = args.b2
        else:
          hap2Bam = outDir+'/HAP2/STARalign/Aligned.out.sorted.bam'
        if args.br:
          refBam = args.br
        else:
          refBam = outDir+'/REF/STARalign/Aligned.out.sorted.bam'
        if not args.c:
          sys.stderr.write("ERROR: rPGA run discover command requires -c parameter \nExample: rPGA run discover -c 1 -v vcf_directory -g annotation.gtf  -o rPGA \n")
          sys.exit()
        chromosome = args.c
        if chromosome.startswith('chr'):
          chromosome = chromosome[3:]

#        if multiprocessing:
#          import multiprocessing
#          CHROMS = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
#          processes = [mp.Process(target=worker, args=(c,)) for c in CHROMS]
#          for p in processes:
#            p.start()
#          for p in processes:
#            p.join()

        else:
          p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting,rnaedit,editFile,gzipped,printall,consensus)
          p.haplotype_specific_junctions()

    elif setting == "alleles":
      if command[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      elif len(command) > 2:
        sys.stderr.write("ERROR: Input arguments are not correct\n")
        sys.exit()
      else:
        discoverJunctions = False        
        writeBam = True
#        if args.v:
#          vcf = args.v
#        else:
#          vcf = open(".rPGAGenotype.yaml").readline().rstrip()
        if not args.v:
          sys.stderr.write("ERROR: rPGA run alleles command requires -v parameter \nExample: rPGA run alleles -c 1 -v vcf_directory -g annotation.gtf  -o rPGA \n")
          sys.exit()
        gtf = ""
        if args.b1:
          hap1Bam = args.b1
        else:
          hap1Bam = outDir+'/HAP1/STARalign/Aligned.out.sorted.bam'
        if args.b2:
          hap2Bam = args.b2
        else:
          hap2Bam = outDir+'/HAP2/STARalign/Aligned.out.sorted.bam'
        refBam = ""
        if not args.c:
          sys.stderr.write("ERROR: rPGA run alleles command requires -c parameter \nExample: rPGA run alleles -c 1 -v vcf_directory -g annotation.gtf  -o rPGA \n")
          sys.exit()
        chromosome = args.c
        if chromosome.startswith('chr'):
          chromosome = chromosome[3:]
        p = DiscoverSpliceJunctions(outDir, vcf, gtf, hap1Bam, hap2Bam, refBam, chromosome, writeBam, discoverJunctions,writeConflicting,rnaedit,editFile,gzipped,printall,consensus)
        p.haplotype_specific_junctions()

    else :
      sys.stderr.write("ERROR: rPGA run -- unknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")
    
