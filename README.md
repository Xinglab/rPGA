##rPGA: RNA-seq Personal Genome-alignment Analyzer
                                          ____  _________
                                    _____/ __ \/ ____/   |
                                   / ___/ /_/ / / __/ /| |
                                  / /  / ____/ /_/ / ___ |
                                 /_/  /_/    \____/_/  |_|

                               ****************************
                               *         V 1.0.1          *
                               ****************************


rPGA is a pipeline to discover  hidden  splicing  variations  by  mapping
personal transcriptomes to personal genomes. In version 1.0.1, rPGA will also
generate allele specific alignments using the "run alleles" option (see Usage).

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

rPGA makes use of python packages pysam and pybedtools to process STAR output.
pysam and pybedtools are their respective dependencies must be installed.


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

Usage: Discover Hidden Splice Junctions
---------------------------------------

Before you can begin working on a project with rPGA, you need to initialize the
project directory. To do this, run:

    $ rPGA init project name

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

rPGA personalize parameters:
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

rPGA mapping parameters:

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

rPGA discover parameters:

     -c CHROM        Chromosome to analyze (required)
     -b              flag to write allele specific bam files
     --conflict      flag to write bam file containing conflicting reads

### Output Files
rPGA run discover outputs 4 files per chromosome:

1. hap1.chrom.specific.bed
2. hap2.chrom.specific.bed
3. hap1hap2.chrom.specific.bed
4. ref.chrom.specific.bed

Columns of each bed file are:

1. chrom
2. junction start
3. junction end
4. name
5. strand
6. splice site usage frequency

The name of each splice junction is in the format J\_R/NC/N3/N5/N35\_SNPid.
- R = junction is in the provided reference annotation
- NC = novel combination of reference 5' and 3' splice sites
- N3 = novel 3'SS and reference 5'SS
- N5 = novel 5'SS and reference 3'SS
- N35 = both 5'SS and 3'SS are novel

SNPid is a comma deliminated list of splice site SNPs.

Usage: Allele Specific Bam Files
--------------------------------

To discover hidden splice junctions, rPGA generates allele specific bam files. If
you are only interested in generating the allele specific bam files, run:

    $ rPGA run alleles

during the last step instead of running the discover function.

rPGA alleles parameters:

     -c CHROM           Chromosome to be analyzed (required)
     --conflict         flag to write bam file containing conflicting reads

Once you have generated allele specific bam files for all 22 autosomal chromosomes
or all 22 autosomes, X, and Y, you can merge them into one allele specific bam file
for each haplotype, hap1.as.bam and hap2.as.bam.

To merge all 22 autosomes run:

    $ rPGA merge auto

To merge all 22 autosome, X, and Y, run:

    $ rPGA merge all



### Generating Allele Specific Bam Files

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
                                           |_______________|
                                                   |
                                                   |
                                      _____________|_____________        
                                     | - Collect reads covering  |
                                     |  heterozygous SNPs        |
                                     | - Perform allele specific |
                                     | 	read assignment          |
                                                   |
                                      _____________|____________
                                 ____|____                  ____|____	 
                                (  Hap 1  )                (  Hap 2  )  
                                ( Specific)                ( Specific)
                                (Alignment)                (Alignment)


To generate allele specific bam files:

1. Personalize reference genome according to a sample's genotype, producing two
personalized genomes, hap1.fa and hap2.fa.

2. Use STAR to align the sample's sequenced reads to the personalized genomes.

3. Collect reads that cover each heterozygous SNP.

4. For each read:
  * Check the position in the read corresponding to the heterozygous SNP.

  * A read is hap1 specific if:
    1. SNP read base matches the hap1 SNP allele
    2. Edit distance to hap1 genome < edit distance to hap2 genome

  * Likewise, a read is hap2 specific if
    1. SNP read base matches the hap2 SNP allele
    2. Edit distance to hap2 genome < edit distance to hap1 genome

  * If a read covers multiple heterozygous SNPs, a majority vote is used. For
      example, if a read covers 3 heterozygous SNPs and 2 match the hap1 allele
      and 1 matches the hap2 allele AND the edit distance to hap1 < edit distance
      to hap2, the read is assigned to hap1.

  * If a read cannot be assigned to either hap1 or hap2 according to the above
      rules, it is considered "conflicting" and is not assigned to either haplotype.
      To output such conflicting reads, use the --conflict option when running
      "discover" or "alleles".

5. Write all hap1 and hap2 specific reads to hap1.as.bam and hap2.as.bam, respectively.   



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
Stein S, Lu Z, Bahrami-Samani E, Park JW, Xing Y. Discover hidden splicing variations by mapping personal transcriptomes to personal genomes. Nucleic Acids Res. 2015 Dec 15;43(22):10612-22.

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
