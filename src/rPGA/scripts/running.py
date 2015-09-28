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

################################################################################
##            PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                ##
################################################################################

<<<<<<< HEAD
def STAR_create_genome(project, genome, gnme):

  # build up options
  opts = ""
  opts += (" --runMode genomeGenerate")  # not tracked; magic number
  opts += (" --genomeDir " + str(project) + "/" + str(gnme) + "/" + "STARindex")
  opts += (" --genomeFastaFiles " + str(genome))  # not tracked; magic number
  opts += (" --runThreadN 8")      # not tracked; magic number

  env_cpy = os.environ.copy()

  commandSTAR = ("STAR" + " " + opts)

  print commandSTAR

  # fire off the sub-process
  processSTAR = subprocess.Popen(commandSTAR, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, env=env_cpy)

  processSTAR.communicate()

  return processSTAR.wait()

def STAR_perform_mapping(project, gnme, seqs):

  # build up options
  opts = ""
  opts += (" --genomeDir " + str(project) + "/" + str(gnme) + "/" + "STARindex")
  opts += (" --readFilesIn " +  str(seqs))
  opts += (" --runThreadN 8 --outFilterMultimapNmax 20")
  opts += (" --outFilterMismatchNmax 0")
  opts += (" --outFileNamePrefix " + str(project))
  opts += ("/" + str(gnme) + "/" + "STARalign ")
  opts += ("--outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonical")
  opts += (" --alignIntronMax 300000 --outSJfilterOverhangMin -1 8 8 8")

  env_cpy = os.environ.copy()

  commandSTAR = ("STAR" + " " + opts)

  print commandSTAR

  # fire off the sub-process
  processSTAR = subprocess.Popen(commandSTAR, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, env=env_cpy)

  processSTAR.communicate()

  return processSTAR.wait()

################################################################################
#             PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                 #
################################################################################

def isHelpString(s):
  norm = s.strip().lower()
  return norm == "help" or norm == "--help"

=======
>>>>>>> parent of 53b0495... Adding the code for running STAR
class PersonalizeGenome :
  def __init__(self, outDir, vcf, ref, hap1Ref, hap2Ref):
    self._outDir = outDir
    self._vcf = vcf
    self._ref = ref
    self._hap1Ref = hap1Ref
    self._hap2Ref = hap2Ref

  def readReference(self):
    # read in hg19 reference file
    chroms = defaultdict(lambda: defaultdict(list)) #chroms[chrom] = [chromosome sequence]
    ref_in = open(self._ref)
    for line in ref_in:
        if line[:1] ==">":
            chroms[line] = []
            key = line
        else:
            line = line.rstrip();
            for base in line:
                chroms[key].append(base)
    ref_in.close()
    return chroms


  def personalizeGenome(self): ## personalize reference genome
    # read in hg19 reference file
    hap1 = self.readReference()

    #loop through vcf file and change reference accordingly
    vcf_in = open(self._vcf)
    for line in vcf_in:
        #skip beginning header lines that begin with "#"
        if line.startswith("#"):
            continue;
        #if SNP is on hap1 (1|0 or 1|1)
        elif "1|" in line:
            line = line.rstrip();
            data = line.split("\t");
            if len(data[3])==1 and len(data[4])==1:
                CHROM = ">chr"+data[0]+"\n";
                POS = int(data[1]);
                #replace REF base with ALT base
                hap1[CHROM][POS-1] = data[4];
    vcf_in.close()
    # print annotated hap1 to file
    FOUT = open(self._hap1Ref,'w')
    for k in hap1.keys():
        FOUT.write(str(k))
        FOUT.write(str("".join(hap1[k])))
        FOUT.write("\n")
    FOUT.close()
    hap1.clear()
    # read in hg19 reference file
    hap2 = self.readReference()
    #loop through vcf file and change reference accordingly
    vcf_in = open(self._vcf)
    for line in vcf_in:
        #skip beginning header lines that begin with "#"
        if line.startswith("#"):
            continue;
        #if SNP is on hap2 (0|1 or 1|1)
        elif "|1" in line:
            line = line.rstrip();
            data = line.split("\t");
            if len(data[3])==1 and len(data[4])==1:
                CHROM = ">chr"+data[0]+"\n";
                POS = int(data[1]);
                #replace REF base with ALT base
                hap2[CHROM][POS-1] = data[4]
    vcf_in.close()
    # print annotated hap2 to file
    FOUT = open(self._hap2Ref,'w')
    for k in hap2.keys():
        FOUT.write(str(k))
        FOUT.write(str("".join(hap2[k])))
        FOUT.write("\n")
    FOUT.close()
    hap2.clear()

### end of personalizing genome ###
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
            "                      ADDING A GENOME                       \n" +\
            "------------------------------------------------------------\n" +\
            "For each genome, rPGA needs to know where to find:          \n" +\
            "A fasta file with the full genome, one chromosome per       \n" +\
            "   sequence.                                                \n" +\
            "To add a genome called hg19, where $ is your prompt:        \n" +\
            "                                                            \n" +\
            "Download the reference sequences for the species from       \n" +\
            "http://hgdownload.cse.ucsc.edu/downloads.html.              \n" +\
            "Concatenate all the files from the different chromosome     \n" +\
            "into one single file. For example:                          \n" +\
            "                                                            \n" +\
            "$ cat *.fa > ~/rPGAGenomes/hg19/hg19.fa                     \n" +\
            "                                                            \n" +\
            "and then run:                                               \n" +\
            "                                                            \n" +\
            "$ rPGA genomes add /path/to/genome                          \n" +\
            "                                                            \n" +\
            "for more information about the add command, where $ is your \n" +\
            "prompt, run:                                                \n" +\
            "                                                            \n" +\
            "$ rPGA genomes add help                                     \n"


  if len(args) == 0 or (len(args) == 1 and isHelpString(args[0])) :
    sys.stderr.write(helpStr + "\n\n")
  else :
    command = args[0].strip().lower()
    if command == "personalize" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 1 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
      else :
        outDir = open(".rPGAProject.yaml").readline().rstrip()
        vcf = open(".rPGAGenotype.yaml").readline().rstrip()
        ref = open(".rPGAGenome.yaml").readline().rstrip()
        hap1Ref = os.path.join(outDir, "hap1.fa")
        hap2Ref = os.path.join(outDir, "hap2.fa")
        p = PersonalizeGenome(outDir, vcf, ref, hap1Ref, hap2Ref)
        p.personalizeGenome()
    elif command == "mapping" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 1 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
#      else :
#        genome = settings(args[1])
    elif command == "discover" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 1 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
#      else :
#        genome = settings(args[1])
    else :
      sys.stderr.write("rPGA genomes -- unnknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")
