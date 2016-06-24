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

# standard python includes
import sys,os,argparse

#rPGA includes
import personalize
import mapping
import discover
import assign
import splicing

def main():
  ## this function performs dispatch on the command line that a user proves. Look at the first argument, which specifies he functio nthe user wants rPGA to perform, then dispatch to the appropriate module.

  parser = argparse.ArgumentParser(description='rPGA')
  parser.add_argument('command',nargs='*') ## command user wants to run


  ## multiple modules
  parser.add_argument('--gz',help='flag denoting gzipped reads',action='store_true') ## gzipped reads (for mapping) or VCF file (for everything else) 
  parser.add_argument('-e',help="file containing RNA editing positions, downloaded from RADAR") #for
  parser.add_argument('--rnaedit',help="flag to check for RNA editing events, must also provide an RNA editing file usng -e parameter",action="store_true")
  parser.add_argument('-v',help="VCF genotype directory") ## or vcf file if just providing one
  parser.add_argument('-g',help="GTF annotation file")
  parser.add_argument('-o',help="output directory")
  parser.add_argument('--nmask',help="flag to N-mask reference genome and align to that",action="store_true")

  ## personalize parameters
  parser.add_argument('-r',help="reference fasta file")

  ## mapping parameters
  parser.add_argument('-N',help="number of mismatches per read pair")
  parser.add_argument('-T',help="num threads for STAR alignment")
  parser.add_argument('-M',help="max number of multiple alignments in STAR mapping")
  parser.add_argument('-s',help="fastq read file(s), comma deliminated if paired end")
  parser.add_argument('-r1',help="hap1 reference genome file")#personal genome, if providing own
  parser.add_argument('-r2',help="hap2 reference genome file") # personal genome, if providing own
  parser.add_argument('-rn',help="N-masked reference genome file") # if providing own
  parser.add_argument('--readlength',help="read length, used to generate genome index (see STAR parameter sjdbOverhang)")
  parser.add_argument('--genomedir',help="STAR genome index directory (if user is providing own genome index to use for mapping)")
  parser.add_argument('--hap',help="only map reads to hap1 and hap2 personal genomes",action="store_true")
  ## discover/assign parameters
  parser.add_argument('-c',help='chromsome') #if not provided, assume run all chromsosomes
  parser.add_argument('--conflict',help='flag to print conflicting reads',action='store_true')
  parser.add_argument('-b1',help="haplotype 1 bam file") ## alignments to personal genomes, if providing own  
  parser.add_argument('-b2',help="haplotype 2 bam file")## alignments to personal genomes, if providing own  
  parser.add_argument('-br',help="reference alignment bam file")
  parser.add_argument('--nomerge',help="flag to not merge chromosome bam files",action="store_true")
  

  ## splicing parameters
  parser.add_argument('--hap1Bam', help="provide own hap1 allele specific bam file")
  parser.add_argument('--hap2Bam', help="provide own hap2 allele specific bam file")
  parser.add_argument('--asdir',help='directory containing AS events')
  parser.add_argument('--MXE')
  parser.add_argument('--SE')
  parser.add_argument('--RI')
  parser.add_argument('--A3SS')
  parser.add_argument('--A5SS')
  parser.add_argument('--samples',help="sample file for merging")
  parser.add_argument('--merge',help="merge count files for multiple samples",action="store_true")
  parser.add_argument('--pos2id',help="flag to change snp pos to snp id in count files",action="store_true")
  parser.add_argument('--haplotype',help="flag to do hap1 vs hap2", action="store_true")
  parser.add_argument('--anchorlength', help="anchor length used for STAR alignment (default=8)")
  parser.add_argument('--mxeDetail',help="merge mxe detailed read counts",action="store_true")

#  parser.add_argument('--exonreads',help="to include exon reads")
  args = parser.parse_args()
  command = args.command

  if ((args.rnaedit and not args.e) or (args.e and not args.rnaedit)):
    sys.stderr.write("rPGA2: if --rnaedit flag is used, you must also provide a file containing RNA editing locations using -e parameter \n")
    sys.exit()
  if len(command)==0:
    sys.stderr.write('rPGA2: need a command - personalize, mapping, discover, alleles, splicing\n')
    sys.exit()
  if (not args.o):
    sys.stderr.write("rPGA2: -o outDir parameter is required\n")
    sys.exit()
  elif command[0]=="personalize":
    personalize.main(args)
  elif command[0]=="mapping":
    mapping.main(args)
  elif command[0]=="discover":
    discover.main(args)
  elif command[0]=="assign":
    assign.main(args)
  elif command[0]=="splicing":
    splicing.main(args)
  else:
    sys.stderr.write("rPGA: I don't recognise the option '" + command[0] +"'\n")
    sys.exit()
