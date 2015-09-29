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
            "                      ADDING A JUNCTION FILE                \n" +\
            "------------------------------------------------------------\n" +\
            "For an individual, rPGA needs to know where to find:        \n" +\
            "A vcf file with the genotype                                \n" +\
            "To add a genotype called hg19, where $ is your prompt:      \n" +\
            "                                                            \n" +\
            "$ rPGA junctions add /path/to/junctions                     \n" +\
            "                                                            \n" +\
            "for more information about the add command, where $ is your \n" +\
            "prompt, run:                                                \n" +\
            "                                                            \n" +\
            "$ rPGA genotype add help                                    \n"


  if len(args) == 0 or (len(args) == 1 and isHelpString(args[0])) :
    sys.stderr.write(helpStr + "\n\n")
  else :
    command = args[0].strip().lower()
    if command == "add" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 2 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
      else :
        dest_fn = open(".rPGAJunctions.yaml", "w")
        dest_fn.write(args[1] + "\n")
        dest_fn.close()
    else :
      sys.stderr.write("rPGA genomes -- unnknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")
