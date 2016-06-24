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

import sys,os,gzip,subprocess,re,logging,time,datetime,commands,argparse,random
from collections import defaultdict

## class definitions

class PersonalizeGenome:
  def __init__(self, outDir, vcf, ref, hap1Ref, hap2Ref,rnaedit,editFile,gzipped,nRef,nmask):
    self.outDir = outDir
    self.vcf = vcf
    self.ref = ref
    self.hap1Ref = hap1Ref
    self.hap2Ref = hap2Ref
    self.rnaedit = rnaedit
    self.edit = editFile
    self.report = os.path.join(outDir,'report.personalize.txt')
    self.gzipped = gzipped
    self.nRef = nRef
    if os.path.isdir(vcf):
      self.vcfdir = True
    else:
      self.vcfdir = False
    self.nmask = nmask 
  
  def read_in_edit(self):
    sys.stdout.write("Reading in edit file: " + self.edit + "\n")
    edit_in = open(self.edit)
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
 
  def read_reference(self):
    sys.stdout.write("Reading in reference " + self.ref + "\n")
    f = defaultdict(list)
    ref_in = open(self.ref)
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

  def read_in_all_vcf_files(self):
    V1,V2 = defaultdict(lambda: defaultdict(list)),defaultdict(lambda: defaultdict(list))
    if self.vcfdir:
      for f in os.listdir(self.vcf):
        vcf_file = os.path.join(self.vcf,f)
        self.read_in_vcf(vcf_file,V1,V2)
    else:
      vcf_file = self.vcf
      self.read_in_vcf(vcf_file,V1,V2)
    return V1,V2

  def read_in_vcf(self,vcf_fn,v1,v2):
    sys.stdout.write('Reading vcf file: ' + vcf_fn +' \n')
    VCF_HEADER = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    if vcf_fn.split('.')[-1]=="gz":
      vcf_in = gzip.open(vcf_fn)
    else:
      vcf_in = open(vcf_fn)
    for line in vcf_in:
      if line.startswith('#'):
        continue
      result = {}
      fields = line.rstrip().split()
      for i,col in enumerate(VCF_HEADER):
        result[col] = fields[i]
      infos = [x for x in result['INFO'].split(';') if x.strip()]
      for i in infos:
        if '=' in i:
          key,value = i.split('=')
          result[key] = value
      geno = result['SAMPLE'].split(':')[0]
      if ('|' in geno):
        g1 = int(geno.split('|')[0])
        g2 = int(geno.split('|')[1])
        alt = result['ALT'].split(',')
        alleles = [result['REF']] + alt
        if all([re.match(r'[ACGT]',i) for i in alleles]) and all([len(i)==1 for i in alleles]):
          if (g1 != 0):
            v1[result['CHROM']][int(result['POS'])] = [alleles[0],alleles[g1]]
          if (g2 != 0):
            v2[result['CHROM']][int(result['POS'])] = [alleles[0],alleles[g2]]
    vcf_in.close()
    return

  def read_in_all_vcf_nmask(self):
    V = defaultdict(list)
    if self.vcfdir:
      for f in os.listdir(self.vcf):
        vcf_file = os.path.join(self.vcf,f)
        self.read_in_vcf_nmask(vcf_file,V)
    else:
      vcf_file = self.vcf
      self.read_in_vcf(vcf_file,V)
    return V

  def read_in_vcf_nmask(self,vcf_fn,v):
    VCF_HEADER = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE']
    for c in CHROMS:
      if self.gzipped:
        vcf_in = gzip.open(vcf_fn)
      else:
        vcf_in = open(vcf_fn)
      sys.stdout.write('Reading vcf file: '+vcf_fn+'\n')
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
          alt = result['ALT'].split(',')
          alleles = [result['REF']] + alt
          geno = result['SAMPLE'].split(':')[0]
          if all([re.match(r'[ACGT]',i) for i in alleles]) and all([len(i)==1 for i in alleles]):
            if ('|' in geno):
              g1 = int(geno.split('|')[0])
              g2 = int(geno.split('|')[1])
              if (g1 != 0) or (g2!=0):
                v[result['CHROM']].append(int(result['POS']))
            elif ('/' in geno):
              try:
                g1 = int(geno.split('/')[0])
                g2 = int(geno.split('/')[1])
                if (g1 != 0) or (g2!=0):
                   v[result['CHROM']].append(int(result['POS']))
              except:
                 continue
    return v  

  
  def personalize_genome(self):
    sys.stdout.write("Personalizing genome\n")
    
    if self.nmask: # n mask all SNP sites
      nref = self.read_reference()
      vcf = self.read_in_all_vcf_nmask()
      report_out = open(self.report,'w')
      nref_out = open(self.nRef,'w')
      snpCounter = 0
      for chrom in vcf:
        for pos in vcf[chrom]:
          nref['chr'+chrom][int(pos)-1] = 'N'
          snpCounter += 1
      if self.rnaedit: 
         editSites = self.read_in_edit()
         for chrom in editSites:
           if not 'chr'+chrom in nref:
             continue
           for pos in editSites[chrom]:
             editCounter += 1
             nref['chr'+chrom][int(pos)-1] = 'N'
      for chrom in nref:
        nref_out.write('>'+chrom+'\n')
        nref_out.write("".join(nref[chrom])+'\n')
      nref.clear()
      nref_out.close()
      report_out.write('# number of N-masked SNPs: '+ str(snpCounter) + '\n')
    
    elif self.rnaedit: # rna editing-- mask sites
      hap1 = self.read_reference()
      vcf1,vcf2 = self.read_in_all_vcf_files()
      editSites = self.read_in_edit()
      report_out = open(self.report,'w')
      hap1_out = open(self.hap1Ref,'w')
      editCounter = 0
      snpCounter = 0
      for chrom in vcf1:
        if 'chr'+chrom not in hap1:
          continue
        for pos in vcf1[chrom]:
          snpCounter += 1
          ref = vcf1[chrom][pos][0]
          alt = vcf1[chrom][pos][1]
          hap1['chr'+chrom][int(pos)-1] = alt
      for chrom in editSites:
        if 'chr'+chrom not in hap1:
          continue
        for pos in editSites[chrom]:
          editCounter += 1
          hap1['chr'+chrom][int(pos)-1] = 'N'
      for chrom in hap1:
        hap1_out.write('>'+chrom+'\n')
        hap1_out.write(''.join(hap1[chrom])+'\n')
      hap1.clear()
      hap1_out.close()

      report_out.write('# number of hap1 SNPs: ' + str(snpCounter) + '\n')
      report_out.write('# number of editing sites: ' + str(editCounter) + '\n')

      editCounter = 0
      snpCounter = 0
      hap2 = self.read_reference()
      hap2_out = open(self.hap2Ref,'w')
      for chrom in vcf2:
        if 'chr'+chrom not in hap2:
          continue
        for pos in vcf2[chrom]:
          snpCounter += 1
          ref = vcf2[chrom][pos][0]
          alt = vcf2[chrom][pos][1]
          hap2['chr'+chrom][int(pos)-1] = alt
      for chrom in editSites:
        if 'chr'+chrom not in hap2:
          continue
        for pos in editSites[chrom]:
          editCounter += 1
          hap2['chr'+chrom][int(pos)-1] = 'N'
      for chrom in hap2:
        hap2_out.write('>'+ chrom+'\n')
        hap2_out.write(''.join(hap2[chrom])+'\n')
      hap2.clear()
      hap2_out.close()
      
      report_out.write('# number of hap2 SNPs: ' + str(snpCounter) + '\n')
      report_out.write('# number of editing sites: ' + str(editCounter) + '\n')
    
    else: # no rna editing
      hap1 = self.read_reference()
      vcf1,vcf2 = self.read_in_all_vcf_files()
      report_out = open(self.report,'w')
      hap1_out = open(self.hap1Ref,'w')
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
      hap2_out = open(self.hap2Ref,'w')
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


## end of class definitions


## begin main function
def main(args): 
  ## Main entry point for personalize module

  helpStr = "--------------------------------------------- \n"+\
            "                  RUNNING PERSONALIZE          \n"+\
            "--------------------------------------------- \n"+\
            " To run personalize, where $ is your prompt:  \n"+\
            " $rPGA personalize [OPTIONS] -r reference -v vcf -o outdir \n" \
            " ** see README for optional parameters        \n" 

  ## store command line parameters
  command = args.command
  if args.rnaedit:
    editFile = args.e
  else:
    editFile = ""
  namsk = args.nmask
  gzipped = args.gz
  rnaedit = args.rnaedit
  outDir = args.o
  if not os.path.exists(outDir):
    os.makedirs(outDir)
    
  ref = args.r
  vcf = args.v
  gtf = args.g
  nmask = args.nmask
  if args.r1:
    hap1Ref = args.r1
  else:
    hap1Ref = os.path.join(outDir, "hap1.fa")
  if args.r2:
      hap2Ref = args.r2
  else:
      hap2Ref = os.path.join(outDir, "hap2.fa")
  if args.rn:
    nRef = args.rn
  else:
    nRef = os.path.join(outDir, "mask.fa")
  
  if len(command)>1:
    if command[1].strip().lower() == "help" :
      sys.stderr.write(helpStr + "\n\n")
      sys.exit()
  if not args.r:
    sys.stderr.write("ERROR: rPGA  personalize command requires -r parameter \nExample: rPGA personalize -v vcf_directory -r reference.fa  -o rPGA \n")
    sys.exit()
  if not args.v:
    sys.stderr.write("ERROR: rPGA  personalize command requires -v parameter \nExample: rPGA personalize -v vcf_directory -r reference.fa  -o rPGA \n")
    sys.exit()
  if not args.o:
    sys.stderr.write("ERROR: rPGA  personalize command requires -o parameter \nExample: rPGA personalize -v vcf -r reference.fa -o rPGA \n")
    sys.exit()
  p = PersonalizeGenome(outDir, vcf, ref, hap1Ref, hap2Ref,rnaedit,editFile,gzipped,nRef,nmask)
  p.personalize_genome()

  return

  
