##rPGA: RNA-seq Personal Genome-alignment Analyzer
                                          ____  _________
                                    _____/ __ \/ ____/   |
                                   / ___/ /_/ / / __/ /| |
                                  / /  / ____/ /_/ / ___ |
                                 /_/  /_/    \____/_/  |_|

                               ****************************
                               *         V 1.3.4          *
                               ****************************


rPGA is a pipeline to discover  hidden  splicing  variations  by  mapping
personal transcriptomes to personal genomes. As of version 1.0.1, rPGA will also
generate allele specific alignments using the "run alleles" option (see Usage).
As of version 1.1.1, rPGA will N-mask known RNA-editing sites when making the 
personal genomes (see --rnaedit flag in Usage).  

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

Table of contents
-----------------
1. Installing rPGA
2. Usage - Discover hidden splice junctions
3. Usage - Allele specific alignment
4. Contacts and bug reports

1. Installation
------------
All of the below installation instructions assume you have write access to the
installation directory for rPGA; if you want to install into a central
location for all users, you may need to prefix these with sudo.

### Unpack ###
To begin the installation, unpack the distribution and CD into the newly created
directory.

    tar -xf rPGA-1.2.4.tar.gz
    cd rPGA-1.2.4

Or, clone from the github repository.

    git clone https://github.com/Xinglab/rPGA.git
    cd rPGA	

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

2. Usage: Discover Hidden Splice Junctions
------------------------------------------

This is the pipeline to discover splice junctions that are hidden when aligning
an individual's transcript reads to the hg19 reference genome. There are three steps 
in the pipeline, 1. personalizing the genome, 2. aligning reads, and 3. discovering
hidden splice junctions.

### Personalize Genome ###

Inputs:
  1. Reference genome (FASTA)
  2. Directory containing VCF files, one per chromosome

Download the reference sequences for the species from
http://hgdownload.cse.ucsc.edu/downloads.html.

Concatenate all the files from the different chromosome into one single file.
For example

    $ cat \*.fa > ~/rPGAGenomes/hg19/hg19.fa

For any individual, rPGA needs to know where to find a vcf file with the
genotype. vcf files with multiple samples, such as the 1000 Genomes vcf files,
must be split it into separate files for each individual. This is easily done
using bcftools, which can be downloaded from
https://github.com/samtools/bcftools/wiki/HOWTOs.

To extract the individual genotype, if necessary:

    $ bcftools view -s sample_name -v snps -p /path/to/ALL.chr.genotypes.vcf > \
    sample.chr.genotype.vcf

Note rPGA requires a separate vcf file for each chromosome. They should be located
in a directory and named 1.vcf, 2.vcf, 3.vcf,... 

At this step rPGA is ready to make the personalized genome. To do this run this
command:

    $ rPGA run personalize -r reference.fa -v genotype_directory -o output_directory

rPGA personalize options:

     -o *            output directory
     -r *            reference genome
     -v *            VCF directory
     --gz            flag denoting VCF files are gzipped 
     --rnaedit       flag to N-mask rna editing sites
     -e              file containing RNA editing sites, can be downloaded from RADAR
                     (http://rnaedit.com/download)

\* Required parameters

** Note if --rnaedit flag is used, RNA editing file must be provided using -e.
rPGA will change each RNA editing site to an "N" in the personal genomes. The number
and locations of RNA editing sites that overlap  SNPs will be reported in 
report.personalize.txt.

Outputs:
  1. output_directory/hap1.fa (hap1 personal genome)
  2. output_directory/hap2.fa (hap2 personal genome)

### RNA-seq Alignment ###

The next step is to map the sequencing data to the personalized genome using
STAR alignment tool. 

Inputs:
  1. Reference genome (FASTA)
  2. Read sequences (FASTQ), single or paired end

Please note that we chose STAR for faster alignment, but
that obviously comes at a cost, and that is the required memory. Make sure you
have enough memory on your computer that is running the STAR (~ 32 GB,
depending on the options of the mapper).

To align reads to personal genomes:

    $ rPGA run mapping -g reference.fa -s reads_1.fastq,[reads_2.fastq]

rPGA mapping options:

     -r *       reference genome (Fasta)
     -s *       read sequences, either single or paired end
     -o *       output directory 
     -T	        number of threads STAR uses, default is 8
     -M         max number of multiple alignments, default is 20
     -N         max number of read mismatches, default is 3
     --gz       flag denoting sequence reads are gzipped

\* Required parameters

Outputs:
  1. output_directory/HAP1/STARalign/Aligned.out.sorted.bam
  2. output_directory/HAP1/STARalign/Aligned.out.sorted.bam.bai
  3. output_directory/HAP2/STARalign/Aligned.out.sorted.bam
  4. output_directory/HAP2/STARalign/Aligned.out.sorted.bam.bai
  5. output_directory/REF/STARalign/Aligned.out.sorted.bam
  6. output_directory/REF/STARalign/Aligned.out.sorted.bam.bai

### Discover hidden splice junctions ###

The final step is to discover novel junctions using the discover function.

Inputs:
  1. Chromosome to be analyzed
  2. Annotation file (GTF)
  3. Genotype directory containing VCF files

Usage:

    $ rPGA run discover -c CHROM -g annotation.gtf -v genotype_directory -o output

rPGA discover options:

     -c *            Chromosome to analyze (required)
                     Note: can be in the form -c 1 or -c chr1 to analyze chrom 1  
     -g *            Annotation file (GTF)
     -v *            Genotype directory (VCF)
     -o *            Output directory
     --writeBam      flag to write allele specific bam files
     --conflict      flag to write bam file containing conflicting reads
     --rnaedit **    flag to consider rna editing sites 
     -e **           file containing RNA editing sites; can be downloaded from
                     RADAR (www.rnaedit.com/download/)
     --gz            flag denoting VCF genotype files are gzipped
     --printall      flag to include non-haplotype specific reads in bam output file
     --consensus     flag to print consensus BAM file (described below)
     -b1 ***         Haplotype 1 alignment file to personal genome (BAM)
     -b2 ***         Haplotype 2 alignment file to personal genome (BAM)
     -br ***         Reference alignment file to reference genome (BAM)

\* Required parameters

** Note: if --rnaedit flag is used, a file containing RNA editing events must be
provided using -e. In this case, rPGA will disregard heterozygous SNPs that overlap
RNA editing sites when assigning mapped reads to haplotypes. Reads that cover RNA editing 
sites will be printed to hap1.chrom.rnaedit.bam and hap2.chrom.rnaedit.bam

*** Use these if you would like to supply your own personal genome mapping alignment files. 
If rPGA run mapping (previous step) is used there is no need to provide rPGA the alignment 
files, rPGA will use the alignments in HAP1/STARalign, HAP2/STARalign, REF/STARalign. 

Outputs (per chromosome): 

Haplotype specific bed files:
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

SNPid is a comma deliminated list of the splice site SNP ids, which match the 
SNP ids in the given VCF file..

If --rnaedit is used, rPGA will also output bam files containing reads 
that overlap RNA editing sites:
  1. hap1.chrom.rnaedit.bam
  2. hap2.chrom.rnaedit.bam

If --writeBam flag is used, rPGA will output allele specific bam files:
  1. hap1.chrom.bam
  2. hap2.chrom.bam

If --consensus flag is used, rPGA will output one consensus bam file:
  1. consensus.chrom.bam

Usage: Allele Specific Bam Files
--------------------------------

To discover hidden splice junctions, rPGA generates allele specific bam files. Follow
this pipeline if you are only interested in generating the allele specific bam files.


### Personalize Genome ###

Inputs:
  1. Reference genome (FASTA)
  2. Directory containing VCF files, one per chromosome

Download the reference sequences for the species from
http://hgdownload.cse.ucsc.edu/downloads.html.

Concatenate all the files from the different chromosome into one single file.
For example

    $ cat \*.fa > ~/rPGAGenomes/hg19/hg19.fa

For any individual, rPGA needs to know where to find a vcf file with the
genotype. vcf files with multiple samples, such as the 1000 Genomes vcf files,
must be split it into separate files for each individual. This is easily done
using bcftools, which can be downloaded from
https://github.com/samtools/bcftools/wiki/HOWTOs.

To extract the individual genotype, if necessary:

    $ bcftools view -s sample_name -v snps -p /path/to/ALL.chr.genotypes.vcf > \
    sample.chr.genotype.vcf

Note rPGA requires a separate vcf file for each chromosome. They should be located
in a directory and named 1.vcf, 2.vcf, 3.vcf,... 

At this step rPGA is ready to make the personalized genome. To do this run this
command:

    $ rPGA run personalize -r reference.fa -v genotype_directory -o output_directory

rPGA personalize options:

     -o *            output directory
     -r *            reference genome
     -v *            VCF directory
     --gz            flag denoting VCF files are gzipped 
     --rnaedit       flag to N-mask rna editing sites
     -e              file containing RNA editing sites, can be downloaded from RADAR
                     (http://rnaedit.com/download)

\* Required parameters

** Note if --rnaedit flag is used, RNA editing file must be provided using -e.
rPGA will change each RNA editing site to an "N" in the personal genomes. The number
and locations of RNA editing sites that overlap  SNPs will be reported in 
report.personalize.txt.

Outputs:
  1. output_directory/hap1.fa (hap1 personal genome)
  2. output_directory/hap2.fa (hap2 personal genome)

### RNA-seq Alignment ###

The next step is to map the sequencing data to the personalized genome using
STAR alignment tool. 

Inputs:
  1. Reference genome (FASTA)
  2. Read sequences (FASTQ), single or paired end

Please note that we chose STAR for faster alignment, but
that obviously comes at a cost, and that is the required memory. Make sure you
have enough memory on your computer that is running the STAR (~ 32 GB,
depending on the options of the mapper).

To align reads to personal genomes:

    $ rPGA run mapping alleles -g reference.fa -s reads_1.fastq,[reads_2.fastq]

rPGA mapping options:

     -r *       reference genome (Fasta)
     -s *       read sequences, either single or paired end
     -o *       output directory 
     -T	        number of threads STAR uses, default is 8
     -M         max number of multiple alignments, default is 20
     -N         max number of read mismatches, default is 3
     --gz       flag denoting sequence reads are gzipped

\* Required parameters

Outputs:
1. output_directory/HAP1/STARalign/Aligned.out.sorted.bam
2. output_directory/HAP1/STARalign/Aligned.out.sorted.bam.bai
3. output_directory/HAP2/STARalign/Aligned.out.sorted.bam
4. output_directory/HAP2/STARalign/Aligned.out.sorted.bam.bai

### Allele Specific Assignment ###

Finally, assign mapped reads to hap1 or hap2 using alleles function. 

Inputs:
  1. Chromosome to analyze
  2. Directory containing one VCF files

Outputs:
  1. hap1.chrom.bam
  2. hap2.chrom.bam
  3. report.chrom.txt

Usage:
	
    $ rPGA run alleles -c CHROM -v genotype_directory -o output_directory

rPGA alleles options:

     -c *           Chromosome to be analyzed (C or -chrC to analyze chrom C)
     -v *           Genotype directory containing VCF files
     -o *           Output directory
     --conflict     flag to write bam file containing conflicting reads
     --rnaedit **   flag to consider rna editing sites 
     -e	**           file containing RNA editing sites; can be downloaded from
                     RADAR (www.rnaedit.com/download/)
     --gz            flag denoting VCF genotype files are gzipped
     --printall      flag to include non-haplotype specific reads in bam output file
     --consensus     flag to print consensus bam file
     -b1 ***         Haplotype 1 alignment file to personal genome (BAM)
     -b2 ***         Haplotype 2 alignment file to personal genome (BAM)
    
\* Required parameter

** Note: if --rnaedit flag is used, a file containing RNA editing events must be 
provided using -e. In this case, rPGA will disregard heterozygous SNPs that overlap
RNA editing sites when assigning mapped	reads to haplotypes.

*** Use these if you would like to supply your own personal genome mapping alignment files. 
If rPGA run mapping (previous step) is used there is no need to provide rPGA the alignment 
files, rPGA will use the alignments in HAP1/STARalign and HAP2/STARalign.

Once you have generated allele specific bam files for all 22 autosomal chromosomes
or all 22 autosomes and X, you can merge them into one allele specific bam file
for each haplotype, hap1.bam and hap2.bam.

To merge all 22 autosomes run:

    $ rPGA merge auto

To merge all 22 autosome and X:

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

5. Write all hap1 and hap2 specific reads to hap1.bam and hap2.bam, respectively.   

Note, rPGA adds up to 4 SAM tags:

1. HT: denotes which haplotype read originates from (1 or 2)
2. SP: comma deliminated list of heterozygous SNP positions the read covers
3. GT: denotes whether the read contains the reference or alternate allele(s) for the
SNP(s) in SP
4. EP: comma deliminated list of RNA editing positions the read covers (if used —-rnaedit)

### Consensus BAM file ###

The consensus BAM file is a single BAM file that contains the best matching reads based on
the haplotype assignment + non heterozygous reads alignments. It consists of reads obtained 
by the following procedure:

1. Reads that uniquely map to one haplotype or the other.
2. Choose read that maps better to one haplotype (due to heterozygous SNP)
3. Randomly choose a read if the read has the same alignment to both haplotypes.
4. Randomly choose a read if the read is conflicting. 
5. Disregard reads that map to different locations in each haplotype, or are multiply mapped.

These comparisons are made for every read that is aligned. For example, if your initial fastq
file contains 100 million reads, there would be 100 million reads in the consensus BAM file, 
minus the reads discarded in step 5. 

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
