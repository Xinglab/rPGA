##rPGA: RNA-seq Personal Genome-alignment Analyzer

                                               ____  _________
                                         _____/ __ \/ ____/   |
                                        / ___/ /_/ / / __/ /| |
                                       / /  / ____/ /_/ / ___ |
                                      /_/  /_/    \____/_/  |_|
                                  *********************************
                                  *            V 1.0.1            *
                                  *********************************


rPGA is a pipeline to discover  hidden  splicing  variations  by  mapping
personal transcriptomes to personal genomes. In version 1.0.1, rPGA will also 
generate allele specific alignments. 

Overview
------------

                     ________________           ____________________
                    (Reference Genome)         (Genotype Information)  
                             |____________________________|
                                  _________|_________
                                 |  Personal Genome  |
                                           |
                                           |
                                       ____|____           _______________
                                      |   STAR  | <------ (Sequenced Reads)
                            _______________|_______________
                       ____|____       ____|____       ____|____
                      (Reference)     (  Hap 1  )     (  Hap 2  )
                      (Alignment)     (Alignment)     (Alignment)
                           |_______________|_______________|
                                           |
                          _________________|___________________
                         | - Discover Personal Splice junctions|       ______________________
                         | - Calculate Usage Frequencies       | <--- (Known Splice Junctions)
                         | - Compare to Known Splice Junctions |
                                           |
                          _________________|_________________
                     ____|____   ____|____   ____|____   ____|____
                    |  Hap 1  | |  Hap 2  | |Hap 1 & 2| |Reference|
                    |Specific | |Specific | |Specific | |Specific |
                    |Junctions| |Junctions| |Junctions| |Junctions|




Requirements
------------

rPGA is intended to be used in a Unix-based environment. It has been tested
on Mac OS and Linux.

A working installation of Python is required.

rPGA makes use of the STAR software package for alignment.
STAR must be installed. By default, rPGA expects them to be somewhere in your
path, but you if they are not you can specify their locations when running the
configure script (see below).

Installation
------------
All of the below installation instructions assume you have write access to the
installation directory for rPGA; if you want to install into a central
location for all users, you may need to prefix these with sudo.

### Unpack ###
To begin the installation, unpack the distribution and CD into the newly created
directory.

    tar -xf rPGA-0.0.1.tar.gz
    cd rPGA-0.0.1


### Configure ###
 Now type

    ./configure

This will configure the installation. rPGA is written partially in Python.
By default, it will use whichever python interpreter is invoked when you type
'Python' on the command line. If you wish to specify a different version of
Python, you can do this when running the ``configure`` script. Additionally,
if any of the required external software packages listed above are not on your
path, you can specify their location when running ``configure``. For example,
to specify the path to python and all of the exeternal software packages,
you would type the following:

    ./configure --with-python=/some/path/to/python \
                --with-STAR=/some/path/to/STAR

You needn't include them all, just the ones that cannot be found in your path.
The configure script will tell you if any requirements cannot be found, in
which case rPGA cannot be installed.

### Install ###
After configuration, rPGA is installed by typing:

    make install

This will install Python modules for whichever version of Python you are using.

User executables and scripts will be placed into /path/to/rPGA/bin and
/path/to/rPGA/scripts respectively. These cannot be moved, as the pipeline
expects them to remain here. We suggest you modify your PATH environment
variable to include these directories, but it is not required. All of the
below instructions assume these directories are in your PATH variable; if not,
replace the names of the scripts/executable with their full path.

Usage
-----

Before you can begin working on a project with rPGA, you need to initialize the
project directory. To do this, run:

    $ rPGA init project_name

For each genome, rPGA needs to know where to find A fasta file with the full
genome, one chromosome per sequence. To add a genome called hg19, where $ is
your prompt:

Download the reference sequences for the species from
http://hgdownload.cse.ucsc.edu/downloads.html.

Concatenate all the files from the different chromosome into one single file.
For example

    $ cat \*.fa > ~/rPGAGenomes/hg19/hg19.fa

and then run:

    $ rPGA genomes add /path/to/genome

For any individual, rPGA needs to know where to find a vcf file with the
genotype. vcf files with multiple samples, such as the 1000 Genomes vcf files,
must be split it into separate files for each individual. This is easily done
using bcftools, which can be downloaded from
https://github.com/samtools/bcftools/wiki/HOWTOs.

To add a genotype, where $ is your prompt:

First extract the individual genotype, if necessary:

    $ bcftools view -s sample_name -v snps -p /path/to/ALL.chr.genotypes.vcf > \
    sample.chr.genotype.vcf

Note rPGA requires a separate vcf file for each chromosome. They should be located
in a directory and named 1.vcf, 2.vcf, 3.vcf,...

and then run:

    $ rPGA genotype add /path/to/genotype_directory

Where 1.vcf, 2.vcf, 3.vcf are located in genotype_directory

At this step rPGA is ready to make the personalized genome. To do this run this
command:

    $ rPGA run personalize

rPGA personalize options:
     -T		 number of threads to use when building STAR genome, default is 8

The next step is to map the sequencing data to the personalized genome using
STAR alignment tool. To do this, first rPGA needs to know where to find the
sequenced reads. To add the sequence files run:

    $ rPGA sequences add /path/to/sequences.fastq

If you have paired end data add the sequences in one line right after each
other, for example:

    $ rPGA sequences add /path/to/sequences_mate_1.fastq /path/to/sequences_mate_2.fastq

Then you can run:

    $ rPGA run mapping

to perform the mapping. Please note that we chose STAR for faster alignment, but
that obviously comes at a cost, and that is the required memory. Make sure you
have enough memory on your computer that is running the STAR (~ 16-28 GB,
depending on the options of the mapper).

rPGA mapping options:

     -T	     number of threads STAR uses, default is 8
     -M      max number of multiple alignments, default is 20
     -N      max number of read mismatches, default is 3
     -g	     use if sequence reads are gzipped

After this step is done, it is time to discover novel junctions, in order to do
this, rPGA needs to know where to find the known splice junction. To add the
know splice junctions, you can run:

    $ rPGA junctions add /path/to/known_splice_junctions

Finally, rPGA is ready to discover novel splice junctions. To do this run:

    $ rPGA run discover

rPGA discover options:

     -c CHROM 		Chromosome to analyze 
     -b			flag to write allele specific bam files 
     --conflict		flag to write bam file containing conflicting reads 

If you are only interested in generating the allele specific bam files, run:

    $ rPGA run alleles

rPGA alleles options:

     -c CHROM 		Chromosome to be analyzed
     --conflict		flag to write bam file containing conflicting reads

   			

Enjoy!

Contacts and bug reports
------------------------
Yi Xing
yxing@ucla.edu

Shayna Stein
sstein93@ucla.edu

Emad Bahrami-Samani
ebs@ucla.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.

Publication
------------
Stein S, Lu Z, Bahrami-Samani E, Park JW, Xing Y. Discover hidden splicing variations by mapping personal transcriptomes to personal genomes. Nucleic Acids Res. In Press.

Copyright and License Information
---------------------------------
Copyright (C) 2015 University of California, Los Angeles (UCLA)
Shayna R. Stein, Emad Bahrami-Samani, Yi Xing

Authors: Shayna R. Stein, Emad Bahrami-Samani, Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
