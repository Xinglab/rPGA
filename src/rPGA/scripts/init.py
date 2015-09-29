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
import sys, os, shutil

################################################################################
##            PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                ##
################################################################################

def isHelpString(s) :
  norm = s.strip().lower()
  return norm == "help" or norm == "--help"

def main(args) :
  """
    Main entry point for this script.
  """

  helpStr = "------------------------------------------------------------\n" +\
            "                    INITIALISE A PROJECT                    \n" +\
            "------------------------------------------------------------\n" +\
            "Before you can begin working on a project with rPGA, you    \n" +\
            "need to initialise the project directory. To do this, first \n" +\
            "move all of your fastq files into the directory (sym-links, \n" +\
            "or hard-links are acceptable), then (assuming $ is your     \n" +\
            "prompt) run this command:                                   \n" +\
            "                                                            \n" +\
            "$ rPGA init project name                                    \n" +\
            "                                                            \n" +\
            "for more information about the project init command, where  \n" +\
            "$ is your prompt, run:                                      \n" +\
            "                                                            \n" +\
            "$ rPGA init --help                                          \n"

  if len(args) == 0 or (len(args) == 1 and isHelpString(args[0])) :
    sys.stderr.write(helpStr + "\n\n")
  else :
    command = args[0].strip().lower()
    if command == "project" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 2 :
        sys.stderr.write("Thename of the project needs to be specified\n")
        sys.exit()
      else :
        if not os.path.exists(args[1]):
          os.makedirs(args[1])
          os.makedirs(os.path.join(args[1], "HG19"))
          os.makedirs(os.path.join(args[1], "HAP1"))
          os.makedirs(os.path.join(args[1], "HAP2"))
          os.makedirs(os.path.join(args[1], "HG19/temp"))
          os.makedirs(os.path.join(args[1], "HAP1/temp"))
          os.makedirs(os.path.join(args[1], "HAP2/temp"))
          os.makedirs(os.path.join(args[1], "HG19/STARindex"))
          os.makedirs(os.path.join(args[1], "HAP1/STARindex"))
          os.makedirs(os.path.join(args[1], "HAP2/STARindex"))
          os.makedirs(os.path.join(args[1], "HG19/STARalign"))
          os.makedirs(os.path.join(args[1], "HAP1/STARalign"))
          os.makedirs(os.path.join(args[1], "HAP2/STARalign"))
          dest_fn = open(".rPGAProject.yaml", "w")
          dest_fn.write(args[1])
          dest_fn.close()
        else :
          sys.stderr.write("Path already exists\n")
          sys.exit()
    else :
       sys.stderr.write("rPGA genomes -- unnknown command: " + command + "\n")
       sys.stderr.write(helpStr + "\n\n")
