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
import sys, os

# rPGA
import init
import genomes
import genotype
import junctions
import seqs
import running
import argparse
#import running2
import merge
def main() :
  """
    This function just performs dispatch on the command line that a user
    provided. Basically, look at the first argument, which specifies the
    function the user wants rPGA to perform, then dispatch to the appropriate
    module.
  """
  parser = argparse.ArgumentParser(description='rPGA')
  parser.add_argument('command',nargs='*')
  parser.add_argument('-N')
  parser.add_argument('-T')
  parser.add_argument('-c',help='chromsome')
  parser.add_argument('-b',help='flag to write allele specific bam files',action='store_true')
  parser.add_argument('-p',help='multiprocessing flag',action='store_true')
  parser.add_argument('-g',help='flag denoting gzipped reads',action='store_true')
  parser.add_argument('--conflict',help='flag to print conflicting reads',action='store_true')
  parser.add_argument('-M',help="max number of multiple alignments in STAR mapping")
  parser.add_argument('-e',help="file containing RNA editing positions, downloaded from RADAR")
  parser.add_argument('--rnaedit',help="flag to check for RNA editing events, must also provide an RNA editing file usng -e parameter",action="store_true")
  args = parser.parse_args()
  command = args.command
  
  if ((args.rnaedit and not args.e) or (args.e and not args.rnaedit)):
    sys.stderr.write("rPGA: if --rnaedit flag is used, you must also provide a file containing RNA editing locations using -e parameter \n")
    sys.exit()
  if len(command)==0:
    sys.stderr.write("rPGA: need a command - init, genomes, genotype, junctions, sequences, or run \n'" )
    print helpStr
  elif command[0]=="init":
    init.main(command[1:])
  elif command[0] == "genomes" :
    genomes.main(command[1:])
  elif command[0] == "genotype" :
    genotype.main(command[1:])
  elif command[0] == "junctions" :
    junctions.main(command[1:])
  elif command[0] == "sequences" :
    seqs.main(command[1:])
  elif command[0] == "run" :
    running.main(args)
  elif command[0] == "merge":
    merge.main(command[1:])
  else :
    sys.stderr.write("rPGA: I don't recognise the option '" + command[0] +\
                     "'.\n")

"""
  @summary: This is a pipeline called rPGA for processing the ...

"""
