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

import sys,os,gzip,re,logging,time,datetime,commands,argparse,random,shutil
from collections import defaultdict
import pysam
from pysam import Samfile

## function definitions
def STAR_create_genome(project, genome, gnme, threads,gtf,readlength):
  # build up options
  sys.stdout.write("Building genome index \n")
  opts = ""
  opts += (" --runMode genomeGenerate")
  opts += (" --genomeDir " + os.path.join(str(project),str(gnme),"STARindex"))
  opts += (" --genomeFastaFiles " + str(genome))
  opts += (" --runThreadN "+str(threads))
  opts += (" --sjdbGTFfile " + str(gtf))
  opts += (" --sjdbOverhang " + str(readlength))

  env_cpy = os.environ.copy()
  commandSTAR = ("STAR" + " " + opts)

  oFile = open(str(project) + "/mapping_commands.sh","w")
  oFile.write("##### Creating Genome for " + str(gnme) + "#####\n" +\
              commandSTAR + "\n#\n")
  oFile.flush()
  oFile.close()
  status,output=commands.getstatusoutput(commandSTAR)

  return

def STAR_perform_mapping(genomedir,project, gnme, seqs, threads,mismatches,gz, multimapped):
  sys.stdout.write("Aligning reads\n")
  # build up options
  opts = ""
  opts += (" --genomeDir " + str(genomedir))
  opts += (" --readFilesIn " + seqs + ' ')
  if gz:
    opts += (" --readFilesCommand gunzip -c")
  opts += (" --runThreadN " + str(threads))
  opts += (" --outFilterMultimapNmax " + str(multimapped))
  opts += (" --alignEndsType EndToEnd")
  opts += (" --outFilterMismatchNmax " + str(mismatches))
  opts += (" --outFileNamePrefix " + os.path.join(str(project), str(gnme),'STARalign/'))
#  opts += ("/" + str(gnme) + "/" + "STARalign/ ")
  opts += (" --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated")
  opts += (" --alignIntronMax 300000 --outSJfilterOverhangMin 8 8 8 8")
  opts += (" --outSAMtype BAM Unsorted") ## output as bam
  opts += (" --outFilterMultimapScoreRange 0") ## only output top scoring reads

  env_cpy = os.environ.copy()
  commandSTAR = ("STAR" + " " + opts)
  print commandSTAR
  oFile = open(str(project) + "/mapping_commands.sh","a")
  oFile.write("##### Mapping reads for " + str(gnme) + "#####\n" +\
              commandSTAR + "\n#\n")
  oFile.flush()
  oFile.close()
  status,output=commands.getstatusoutput(commandSTAR)

  return

def sam_to_sorted_bam(sam_prefix):
  print "sort bam: ", sam_prefix 
  bam_fn = sam_prefix + '.bam'
  pysam.sort(bam_fn, sam_prefix+'.sorted')
  pysam.index(sam_prefix+'.sorted.bam')
  return

## done defining functions

## begin main

def main(args):
  ## main function for mapping reads to personal genomes
  
  helpStr = "Help!\n"
  if not args.o:
    sys.stderr.write('rPGA2 ERROR: must provide output directory \n\n')
    sys.exit()
  if not args.s:
    sys.stderr.write('rPGA ERROR: must provide read sequence files \n\n')
    sys.exit()
#  if not args.g:
#    sys.stderr.write('rPGA ERROR: must provide gtf file\n')
  

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

  gzipped = args.gz
  outDir = args.o
  nmask = args.nmask
  if not os.path.exists(outDir):
    os.makedirs(outDir)
  ref = args.r
  vcf = args.v
  gtf = args.g
  hap = args.hap
  if args.r1:
    hap1Ref = args.r1
  else:
    hap1Ref = os.path.join(outDir, "hap1.fa")
  if args.r2:
    hap2Ref = args.r2
  else:
    hap2Ref = os.path.join(outDir, "hap2.fa")
  nRef = os.path.join(outDir, "mask.fa")  
  if args.readlength:
    readlength = int(args.readlength)-1
  else:
    readlength = 99
  if len(command)>1:
    sys.stderr.write(helpStr + "\n\n")
    sys.exit()

  seqs = ' '.join((args.s).split(','))
  if len((args.s).split(','))==0 or len((args.s).split(','))>2:
    sys.stderr.write("ERROR: Sequence parameter -s input is  not correct\n Example: rPGA run mappng alleles -s reads_1.fq,reads_2.fq -o rPGA\n")
    sys.exit()
  
  
  if nmask:
    if not os.path.exists(os.path.join(outDir,"MASK/STARindex")):
      os.makedirs(os.path.join(outDir, "MASK/STARindex"))
    if not os.path.exists(os.path.join(outDir,"MASK/STARalign")):
      os.makedirs(os.path.join(outDir, "MASK/STARalign"))
  ##create genome index
    STAR_create_genome(outDir, nRef, "MASK",threads,gtf,readlength)
    genomeDir = outDir + '/MASK/STARindex'
  # map reads
    STAR_perform_mapping(genomeDir,outDir, "MASK", seqs,threads,mismatches,gzipped,multimapped)   
  ## sort bam file
    sam_to_sorted_bam(os.path.join(outDir,'MASK/STARalign/Aligned.out'))
  ## remove unnecessary files
    os.remove(os.path.join(outDir,'MASK/STARalign/Aligned.out.bam'))
    os.remove(os.path.join(outDir,'MASK/STARalign/SJ.out.tab'))
    shutil.rmtree(os.path.join(outDir, "MASK/STARindex"))

  elif hap: # map to hap and hap2 personal genomes
    if not os.path.exists(os.path.join(outDir, "HAP1/STARindex")):
      os.makedirs(os.path.join(outDir, "HAP1/STARindex"))
    if not os.path.exists(os.path.join(outDir, "HAP2/STARindex")):
      os.makedirs(os.path.join(outDir, "HAP2/STARindex"))
    if not os.path.exists(os.path.join(outDir, "HAP1/STARalign")):
      os.makedirs(os.path.join(outDir, "HAP1/STARalign"))
    if not os.path.exists(os.path.join(outDir, "HAP2/STARalign")):
      os.makedirs(os.path.join(outDir, "HAP2/STARalign"))
    if not args.genomedir:
      if not args.g:
        sys.stderr.write('rPGA ERROR: must provide gtf file\n') 
        sys.exit()
      STAR_create_genome(outDir, hap1Ref, "HAP1",threads,gtf,readlength)
      STAR_create_genome(outDir, hap2Ref, "HAP2",threads,gtf,readlength)
      genomeDir1 = os.path.join(outDir, 'HAP1/STARindex')
      genomeDir2 = os.path.join(outDir, 'HAP2/STARindex')
    else:
      genomeDir1, genomeDir2 = (args.genomedir).split(',')
    STAR_perform_mapping(genomeDir1, outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
    STAR_perform_mapping(genomeDir2, outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
    sam_to_sorted_bam(os.path.join(outDir,'HAP1/STARalign/Aligned.out'))
    sam_to_sorted_bam(os.path.join(outDir,'HAP2/STARalign/Aligned.out'))
    os.remove(os.path.join(outDir,'HAP1/STARalign/Aligned.out.bam'))
    os.remove(os.path.join(outDir,'HAP2/STARalign/Aligned.out.bam'))
    os.remove(os.path.join(outDir,'HAP1/STARalign/SJ.out.tab'))
    os.remove(os.path.join(outDir,'HAP2/STARalign/SJ.out.tab'))
    shutil.rmtree(os.path.join(outDir, "HAP1/STARindex"))
    shutil.rmtree(os.path.join(outDir, "HAP2/STARindex"))
  else:
    if not args.r:
      sys.stderr.write("ERROR: rPGA run mapping command requires -r parameter \nExample: rPGA run mapping -r reference.fa -s reads_1.fq,reads_.fq -o rPGA \n")
      sys.exit()
    
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
    if not args.genomedir:
      STAR_create_genome(outDir, ref, "REF",threads,gtf,readlength)
      STAR_create_genome(outDir, hap1Ref, "HAP1",threads,gtf,readlength)
      STAR_create_genome(outDir, hap2Ref, "HAP2",threads,gtf,readlength)

    print "perform STAR mapping"
    if args.genomedir:
      genomeDir1, genomeDir2, genomeDirR = (args.genomedir).split(',')
    else:
      genomeDir1 = os.path.join(outDir, 'HAP1/STARindex')
      genomeDir2 = os.path.join(outDir, 'HAP2/STARindex')
      genomeDirR = os.path.join(outDir, 'REF/STARindex')
    if not args.g:
      sys.stderr.write('rPGA ERROR: must provide gtf file\n')
      sys.exit()
    STAR_create_genome(outDir, hap1Ref, "HAP1",threads,gtf,readlength)
    STAR_create_genome(outDir, hap2Ref, "HAP2",threads,gtf,readlength)
    STAR_create_genome(outDir, ref, "REF",threads,gtf,readlength)
    STAR_perform_mapping(genomeDir1, outDir, "HAP1", seqs,threads,mismatches,gzipped,multimapped)
    STAR_perform_mapping(genomeDir2, outDir, "HAP2", seqs,threads,mismatches,gzipped,multimapped)
    STAR_perform_mapping(genomeDirR, outDir, "REF", seqs,threads,mismatches,gzipped,multimapped)
    sam_to_sorted_bam(os.path.join(outDir,'HAP1/STARalign/Aligned.out'))
    sam_to_sorted_bam(os.path.join(outDir,'HAP2/STARalign/Aligned.out'))
    sam_to_sorted_bam(os.path.join(outDir,'REF/STARalign/Aligned.out'))
#    os.remove(os.path.join(outDir,'HAP1/STARalign/Aligned.out.bam'))
#    os.remove(os.path.join(outDir,'HAP2/STARalign/Aligned.out.bam'))
#    os.remove(os.path.join(outDir,'REF/STARalign/Aligned.out.bam'))
#    os.remove(os.path.join(outDir,'HAP1/STARalign/SJ.out.tab'))
#    os.remove(os.path.join(outDir,'HAP2/STARalign/SJ.out.tab'))
#    os.remove(os.path.join(outDir,'REF/STARalign/SJ.out.tab'))
#    shutil.rmtree(os.path.join(outDir, "HAP1/STARindex"))
#    shutil.rmtree(os.path.join(outDir, "HAP2/STARindex"))
#    shutil.rmtree(os.path.join(outDir, "REF/STARindex"))
    
