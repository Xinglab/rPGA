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
import pysam
from pysam import Samfile
import argparse
################################################################################
##            PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                ##
################################################################################

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
            "                      MERGING BAM FILES                     \n" +\
            "------------------------------------------------------------\n" +\
            "After generating allele specific bam files for each         \n" +\
            "chromosome, rPGA can merge them.                            \n" +\
            "                                                            \n" +\
            "To merge bam files for all 22 autosomal chromosomes:        \n" +\
            "                                                            \n" +\
            "$ rPGA merge auto                                           \n" +\
            "                                                            \n" +\
            "To merge bam files for all 22 autosomes and X and Y:        \n" +\
            "                                                            \n" +\
            "$ rPGA merge all                                            \n"


if len(args) < 1:
    sys.stderr.write(helpStr + "\n\n")
    sys.exit()
else :
    command = args[0].strip().lower()
    if command == "all" :
      if args[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      else :
        outDir = open(".rPGAProject.yaml").readline().rstrip()
        pysam.merge('-f',outDir+'/hap1.as.bam',outDir+'/hap1.1.as.bam',outDir+'/hap1.2.as.bam',outDir+'/hap1.3.as.bam',outDir+'/hap1.4.as.bam',outDir+'/hap1.5.as.bam',outDir+'/hap1.6.as.bam',
                    outDir+'/hap1.7.as.bam',outDir+'/hap1.8.as.bam',outDir+'/hap1.9.as.bam',outDir+'/hap1.10.as.bam',outDir+'/hap1.11.as.bam',outDir+'/hap1.12.as.bam',outDir+'/hap1.13.as.bam',
                    outDir+'/hap1.14.as.bam',outDir+'/hap1.15.as.bam',outDir+'/hap1.16.as.bam',outDir+'/hap1.17.as.bam',outDir+'/hap1.18.as.bam',outDir+'/hap1.19.as.bam',outDir+'/hap1.20.as.bam',
                    outDir+'/hap1.21.as.bam',outDir+'/hap1.22.as.bam',outDir+'/hap1.X.as.bam',outDir+'/hap1.Y.as.bam')
        pysam.merge('-f',outDir+'/hap2.as.bam',outDir+'/hap2.1.as.bam',outDir+'/hap2.2.as.bam',outDir+'/hap2.3.as.bam',outDir+'/hap2.4.as.bam',outDir+'/hap2.5.as.bam',outDir+'/hap2.6.as.bam',
                    outDir+'/hap2.7.as.bam',outDir+'/hap2.8.as.bam',outDir+'/hap2.9.as.bam',outDir+'/hap2.10.as.bam',outDir+'/hap2.11.as.bam',outDir+'/hap2.12.as.bam',outDir+'/hap2.13.as.bam',
                    outDir+'/hap2.14.as.bam',outDir+'/hap2.15.as.bam',outDir+'/hap2.16.as.bam',outDir+'/hap2.17.as.bam',outDir+'/hap2.18.as.bam',outDir+'/hap2.19.as.bam',outDir+'/hap2.20.as.bam',
                    outDir+'/hap2.21.as.bam',outDir+'/hap2.22.as.bam',outDir+'/hap2.X.as.bam',outDir+'/hap2.Y.as.bam')
    elif command == "auto":
      if args[1].strip().lower() == "help" :
        sys.stderr.write(helpStr + "\n\n")
        sys.exit()
      else:
        pysam.merge('-f',outDir+'/hap1.as.bam',outDir+'/hap1.1.as.bam',outDir+'/hap1.2.as.bam',outDir+'/hap1.3.as.bam',outDir+'/hap1.4.as.bam',outDir+'/hap1.5.as.bam',outDir+'/hap1.6.as.bam',
                    outDir+'/hap1.7.as.bam',outDir+'/hap1.8.as.bam',outDir+'/hap1.9.as.bam',outDir+'/hap1.10.as.bam',outDir+'/hap1.11.as.bam',outDir+'/hap1.12.as.bam',outDir+'/hap1.13.as.bam',
                    outDir+'/hap1.14.as.bam',outDir+'/hap1.15.as.bam',outDir+'/hap1.16.as.bam',outDir+'/hap1.17.as.bam',outDir+'/hap1.18.as.bam',outDir+'/hap1.19.as.bam',outDir+'/hap1.20.as.bam',
                    outDir+'/hap1.21.as.bam',outDir+'/hap1.22.as.bam')
        pysam.merge('-f',outDir+'/hap2.as.bam',outDir+'/hap2.1.as.bam',outDir+'/hap2.2.as.bam',outDir+'/hap2.3.as.bam',outDir+'/hap2.4.as.bam',outDir+'/hap2.5.as.bam',outDir+'/hap2.6.as.bam',
                    outDir+'/hap2.7.as.bam',outDir+'/hap2.8.as.bam',outDir+'/hap2.9.as.bam',outDir+'/hap2.10.as.bam',outDir+'/hap2.11.as.bam',outDir+'/hap2.12.as.bam',outDir+'/hap2.13.as.bam',
                    outDir+'/hap2.14.as.bam',outDir+'/hap2.15.as.bam',outDir+'/hap2.16.as.bam',outDir+'/hap2.17.as.bam',outDir+'/hap2.18.as.bam',outDir+'/hap2.19.as.bam',outDir+'/hap2.20.as.bam',
                    outDir+'/hap2.21.as.bam',outDir+'/hap2.22.as.bam')
    else :
      sys.stderr.write("rPGA genomes -- unnknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")