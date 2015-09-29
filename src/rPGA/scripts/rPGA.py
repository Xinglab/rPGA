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
import seqs
import running

def main() :
  """
    This function just performs dispatch on the command line that a user
    provided. Basically, look at the first argument, which specifies the
    function the user wants rPGA to perform, then dispatch to the appropriate
    module.
  """
  if sys.argv[1] == "init" :
    init.main(sys.argv[2:])
  elif sys.argv[1] == "genomes" :
    genomes.main(sys.argv[2:])
  elif sys.argv[1] == "genotype" :
    genotype.main(sys.argv[2:])
  elif sys.argv[1] == "junctions" :
    junctions.main(sys.argv[2:])
  elif sys.argv[1] == "sequences" :
    seqs.main(sys.argv[2:])
  elif sys.argv[1] == "run" :
    running.main(sys.argv[2:])
  else :
    sys.stderr.write("rPGA: I don't recognise the option '" + sys.argv[1] +\
                     "'.\n")

"""
  @summary: This is a pipeline called rPGA for processing the ...

"""
