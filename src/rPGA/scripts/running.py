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



################################################################################
#             PRELIMINARY COMMAND LINE PROCESSING AND DISPATCH                 #
################################################################################

class PersonalizeGenome :
  def __init__(self, outDir, vcf, gtf, ref, hap1Ref, hap2Ref):
    self._outDir = outDir
    self._vcf = vcf
    self._gtf = gtf
    self._ref = ref
    self._hap1Ref = hap1Ref
    self._hap2Ref = hap2Ref

  def readReference(self):
    # read in hg19 reference file
    chroms = defaultdict(lambda: defaultdict(list))
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

def STAR_create_genome(self, gnme):

  # build up options
    opts = ""
    opts += (" --runMode genomeGenerate")  # not tracked; magic number
    opts += (" --genomeDir " + str(self._outDir) + "/" + str(gnme) + "/" + "STARindex")
    opts += (" --genomeFastaFiles " + str(self._ref))  # not tracked; magic number
    opts += (" --runThreadN 8")      # not tracked; magic number

    env_cpy = os.environ.copy()
    commandSTAR = ("STAR" + " " + opts)

    oFile = open(str(self._outDir) + "/mapping_commands.sh","w")
    oFile.write("##### Creating Genome for " + str(gnme) + "#####\n" +\
                commandSTAR + "\n#\n")
    oFile.flush()
    oFILE.close()
    status,output=commands.getstatusoutput(commandSTAR)

    return 1

#    fire off the sub-process
#    processSTAR = subprocess.Popen(commandSTAR, shell=True,
#                                   stdout=subprocess.PIPE,
#                                   stderr=subprocess.PIPE, env=env_cpy)

#  processSTAR.communicate()

#  return processSTAR.wait()

  def STAR_perform_mapping(self, gnme, seqs):

    # build up options
    opts = ""
    opts += (" --genomeDir " + str(self._outDir) + "/" + str(gnme) + "/" + "STARindex")
    opts += (" --readFilesIn " +  str(seqs))
    opts += (" --runThreadN 8 --outFilterMultimapNmax 20")
    opts += (" --outFilterMismatchNmax 0")
    opts += (" --outFileNamePrefix " + str(self._outDir))
    opts += ("/" + str(gnme) + "/" + "STARalign ")
    opts += ("--outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonical")
    opts += (" --alignIntronMax 300000 --outSJfilterOverhangMin -1 8 8 8")

    env_cpy = os.environ.copy()
    commandSTAR = ("STAR" + " " + opts)
    oFile = open(str(self._outDir) + "/mapping_commands.sh","w")
    oFile.write("##### Creating Genome for " + str(gnme) + "#####\n" +\
                commandSTAR + "\n#\n")
    oFile.flush()
    oFILE.close()
    status,output=commands.getstatusoutput(commandSTAR)

    return 1

#    fire off the sub-process
#    processSTAR = subprocess.Popen(commandSTAR, shell=True,
#                                   stdout=subprocess.PIPE,
#                                   stderr=subprocess.PIPE, env=env_cpy)

#    processSTAR.communicate()

#    return processSTAR.wait()


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

  def listToString(x):
    rVal = '';
    for a in x:
      rVal += a + ' ';
    return rVal;

  def haplotypeSpecificSam(self): ## extract haplotype specific sam files
    ## extract junction reads for hg19, hap1, and hap2
    oFile = open(str(self._outDir) + "/commands.sh","a")
    hgPath = str(self._outDir) + "/HG19"
    h1Path = str(self._outDir) + "/HAP1"
    h2Path = str(self._outDir) + "/HAP2"
    hgTemp = str(self._outDir) + "/HG19/temp"
    h1Temp = str(self._outDir) + "/HAP1/temp"
    h2Temp = str(self._outDir) + "/HAP2/temp"

    for i in ['HG19','HAP1','HAP2']:
        cmd = 'awk \'$6~"S" {print $1}\' ' + self._outDir + '/' + i + '/STARalign/Aligned.out.sam > '+self._outDir+'/'+i+'/temp/truncatedIDs.txt; '
        cmd += 'awk \'NR==FNR{a[$1]++;next}; a[$1]<=0 {print $1} \' ' + self._outDir+'/'+i+'/temp/truncatedIDs.txt '+self._outDir+'/'+i+'/STARalign/Aligned.out.sam  > ' + self._outDir+'/'+i+'/temp/mappedIDs.txt'
        oFile.write('##### getting junction reads for ' + i + '#####\n'+cmd+'\n#\n')
        oFile.flush()
        status,output=commands.getstatusoutput(cmd)


    ## extract hap1 specific junction reads
    cmd = 'awk \' NR==FNR{a[$1]++;next}; a[$1]<=0 \' '+hgTemp+'/mappedIDs.txt '+h1Temp+'/mappedIDs.txt > '+h1Temp+'/mappedIDs.specific.txt ;'
    cmd += 'awk \'NR==FNR{a[$1]++;next}; (a[$1]>0 && $6~"N") \' '+h1Temp +'/mappedIDs.txt '+h1Path+'/STARalign/Aligned.out.sam > '+h1Temp+'/junction.sam ; '
    cmd += 'awk \'NR==FNR{a[$1]++;next}; (a[$1]>0 && $6~"N") \' '+h1Temp +'/mappedIDs.specific.txt '+h1Path+'/STARalign/Aligned.out.sam > '+h1Temp+'/junction.specific.sam '
    oFile.write('### getting hap1 specific junction reads ###\n' + cmd + '\n#\n')
    oFile.flush()
    status,output=commands.getstatusoutput(cmd)

    ## extract hap2 specific junction reads
    cmd = 'awk \' NR==FNR{a[$1]++;next}; a[$1]<=0 \' '+hgTemp+'/mappedIDs.txt '+h2Temp+'/mappedIDs.txt > '+h2Temp+'/mappedIDs.specific.txt ;'
    cmd += 'awk \'NR==FNR{a[$1]++;next}; (a[$1]>0 && $6~"N") \' '+h2Temp +'/mappedIDs.txt '+h2Path+'/STARalign/Aligned.out.sam > '+h2Temp+'/junction.sam; '
    cmd += 'awk \'NR==FNR{a[$1]++;next}; (a[$1]>0 && $6~"N") \' '+h2Temp +'/mappedIDs.specific.txt '+h2Path+'/STARalign/Aligned.out.sam > '+h2Temp+'/junction.specific.sam '
    oFile.write('### getting hap2 specific junction reads ###\n' + cmd + '\n#\n')
    oFile.flush()
    status,output=commands.getstatusoutput(cmd)

    ## extract hg19 specific junction reads
    cmd = 'awk \' NR==FNR{a[$1]++;next}; a[$1]<=0 \' '+h1Temp+'/mappedIDs.txt '+hgTemp+'/mappedIDs.txt > '+hgTemp+'/mappedIDs.specific.temp.txt ;'
    cmd += 'awk \' NR==FNR{a[$1]++;next}; a[$1]<=0 \' '+h2Temp+'/mappedIDs.txt '+hgTemp+'/mappedIDs.specific.temp.txt > '+hgTemp+'/mappedIDs.specific.txt ; '
    cmd += 'awk \'NR==FNR{a[$1]++;next}; (a[$1]>0 && $6~"N") \' '+hgTemp +'/mappedIDs.txt '+hgPath+'/STARalign/Aligned.out.sam > '+hgTemp+'/junction.sam ; '
    cmd += 'awk \'NR==FNR{a[$1]++;next}; (a[$1]>0 && $6~"N") \' '+hgTemp +'/mappedIDs.specific.txt '+hgPath+'/STARalign/Aligned.out.sam > '+hgTemp+'/junction.specific.sam '
    oFile.write('### getting hap2 specific junction reads ###\n' + cmd + '\n#\n')
    oFile.flush()
    status,output=commands.getstatusoutput(cmd)

    return
## end of haplotypeSpecificSam ##

  def junctionCount(slef,sam_fn):
    ## count number of total and haplotype specific distinct reads that span splice junctions and add to SJ.out files
    readsMapped = defaultdict(list) # key=junction, value = list of reads mapping to junction(by start pos)
    nnR = defaultdict(int) #key = junction, value = number of distinct reads mapped to junction

    for line in open(sam_fn):
        line = line.rstrip()
        if ((not line.startswith('@')) and ('N' in line.split()[5]) and ('S' not in line.split()[5])):
            # not a header line, and read contains a junction and no truncation
            ele = line.split()
            cigar = ele[5] # format xMyNzM
            exon_lengths = re.findall(r"(\d+M)",cigar)
            intron_lengths = re.findall(r"(\d+N)",cigar)
            pos = int(ele[3]) #read start location
            # iterate through exons/introns to extract junction coordinates that read spans
            for index,item in enumerate(intron_lengths):
                jS = pos + int(exon_lengths[index][:-1])
                jE = jS + int(intron_lengths[index][:-1]) - 1
                key = str(ele[2]) + '_' + str(jS) + '_' + str(jE) #chr:start:end
                readsMapped[key].append(ele[3]) #append read start location
                pos = jE
    # end of looping through sam file

    #determine nnR
    for k in readsMapped.keys():
        num = len(set(readsMapped[k]))
        nnR[k] = num
    return nnR

## end of junctionCount ##

  def distinctSJOut(self):
    # count number of distinct reads that align to each junction
    # count the total distinct and haplotype specific distinct
    for i in ['HG19','HAP1','HAP2']:
        total = self.junctionCount(self._outDir+'/'+i+'/temp/junction.sam') ## all junction read counts
        specific = self.junctionCount(self._outDir+'/'+i+'/temp/junction.specific.sam') ## distinct junction read counts
        b_fn = open(self._outDir+'/'+i+'/temp/junctions.distinct.bed','w')
        for line in open(self._outDir +'/'+i+'/STARalign/SJ.out.tab'):
            line = line.rstrip()
            ele = line.split()
            chrom = ele[0]
            jS = ele[1]
            jE = ele[2]
            key = chrom + '_' + str(jS) + '_' + str(jE)
            b_fn.write(line + '\t' + str(total[key]) + '\t' + str(specific[key]) + '\n')
        b_fn.close()

## end of distinctSJOut ##

  def readVCF(self):
    vcf_dict = defaultdict(lambda: defaultdict(str)) # vcf[chrom][pos] = snpid
    for line in open(self._vcf):
        if line.startswith('#'):
            continue # skip header lines
        elif ('|' in line): #snp is phased
            ele = line.split()
            if 'chr' in ele[0]:
                chrom = ele[0]
            else:
                chrom = 'chr' + str(ele[0])
            pos = int(ele[1])
            snpid = ele[2]
            if (ele[4] in ['A','C','G','T']): # alt base is determined
                vcf_dict[chrom][pos] = snpid
    return vcf_dict

  def hapSpecific(self,bed_fn,v):
    ## store junctions with >= 2 haplotype specific reads AND a snp in the splice site
    d = defaultdict(list) #d[chr:s:e] = [snpid(s)]
    for line in open(bed_fn):
        line = line.rstrip()
        ele = line.split()
        if (int(ele[10])>=2): #at least 2 haplotype specific reads span junction
            chrom = ele[0]
            jS = int(ele[1])
            jE = int(ele[2])
            strand = ele[3]
            k = chrom + '_' + str(jS) + '_' + str(jE) + '_' + strand
            for num in range(0,2): #check for SNPs in 5'SS
                if ((jS + num) in v[chrom]):
                    d[k].append(v[chrom][jS+num])
            for num in range(-1,1): # check for SNPs in 3'SS
                if ((jE + num) in v[chrom]):
                    d[k].append(v[chrom][jE+num])
    return d

  def readGTF(self):
    exon_start = defaultdict(list) #stored exon_start positions. key = chrom, value = list of 3' positions
    exon_end = defaultdict(list) #stores exon end potitions. key = chrom, value = list of 5' positions

    for line in open(self._gtf): #for each line
        if line.startswith('#'): #comment line
            continue
        else:
            if 'exon' in line.split()[2]: #entry is information about an exon
                chrom = line.split()[0] #chromosome
                start = int(line.split()[3]) #exon start pos
                end = int(line.split()[4]) #exon end pos
                exon_start[chrom].append(start)
                exon_end[chrom].append(end)
    return exon_start,exon_end


  def calculateFrequency(self,specific, fn,gS,gE):
    ## calculate specific junction frequencies
    ## specific is list of haplotype specific junctions (chr_start_end_strand)
    ## fn is SJ.out.distinct.tab file
    ## gS is exon start positions, gE is exon end positions in gtf file
    f = defaultdict(float) # f[junction] = frequency
    novel = defaultdict(bool)
    ## read in SJ.out file ##
    nR = defaultdict(lambda: defaultdict(dict)) # nR[chrom][start][end] = num_reads
    for line in open(fn):
        line = line.rstrip()
        ele = line.split()
        chrom = ele[0]
        jS = ele[1]
        jE = ele[2]
        n = int(ele[9])
        nR[chrom][jS][jE] = n

    for s in specific:
        j = s.split('_')
        chrom = j[0]
        jS = j[1]
        jE = j[2]
        strand = j[3]
        NR=nR[chrom][jS][jE]
        overlapping_junctions = []
        overlapping_counts = []
        novel[s] = True

        for x in nR[chrom]:
            # x is junction start position
            if ((int(x)-1) in gE[chrom]): # x in gtf file
                if (int(x) < int(jS)):
                    for y in nR[chrom][x]:
                         if ((int(y)+1) in gS[chrom]): # y is in gtf file
                             if (int(y) > int(jS)):
                                 overlapping_junctions.append(chrom+'_'+x+'_'+y)
                                 overlapping_counts.append(int(nR[chrom][x][y]))
                elif (int(x)==int(jS)):
                    for y in nR[chrom][x]:
                        if ((int(y)==int(jE)) and ((int(y)+1) in gS[chrom])):
                            novel[s] = False
                        else:
                            if ((int(y)+1) in gS[chrom]):
                                overlapping_junctions.append(chrom+'_'+x+'_'+y)
                                overlapping_counts.append(int(nR[chrom][x][y]))
                elif (int(x)<int(jE)):
                    for y in nR[chrom][x]:
                        if ((int(y)+1) in gS[chrom]):
                            overlapping_junctions.append(chrom+'_'+x+'_'+y)
                            overlapping_counts.append(int(nR[chrom][x][y]))

        freq = NR/float(NR + sum(overlapping_counts))
        f[s] = freq
    return f, overlapping_junctions,overlapping_counts,novel

  def printBed(self,specific,snpid,freq,novel,o_fn):
    ## print output ##
    OUT = open(o_fn,'w')
    i = 1
    OUT.write('#chrom\tstart\tend\tname_SNPid_novel\tfrequency\tstrand\n')
    for j in sorted(specific):
        chrom = j.split('_')[0]
        jS = j.split('_')[1]
        jE = j.split('_')[2]
        strand = j.split('_')[3]
        if strand=="1":
            strand = '+'
        else:
            strand = '-'
        s = ','.join(snpid[j])
        f = freq[j]
        n = novel[j]
        name = 'J'+str(i)+'_'+s+'_'+str(n)

        OUT.write(chrom+'\t'+jS+'\t'+jE+'\t'+name+'\t'+str(f)+'\t'+strand+'\n')

    return
## end of printBed ##

  def haplotypeSpecificJunctions(self):
    ## extract haplotype specific junctions (hap1, hap2, hap1hap2, hg19)
    ## junction files ##
    h1_fn = os.path.join(self._outDir, "HAP1/temp/junctions.distinct.bed")
    h2_fn = os.path.join(self._outDir, "HAP2/temp/junctions.distinct.bed")
    hg_fn = os.path.join(self._outDir, "HG19/temp/junctions.distinct.bed")
    VCF = self.readVCF()
    ## store junctions with >= 2 haplotype specific reads AND a snp in the splice site
    h1_j = self.hapSpecific(h1_fn,VCF)
    h2_j = self.hapSpecific(h2_fn,VCF)
    hg_j = self.hapSpecific(hg_fn,VCF)

    ## get haplotype specific junctions ##
    hap1_specific = [ x for x in h1_j if ((x not in h2_j) and (x not in hg_j))]
    hap2_specific = [ x for x in h2_j if ((x not in h1_j) and (x not in hg_j))]
    hap1hap2_specific = [ x for x in h1_j if ((x in h2_j) and (x not in hg_j))]
    hg19_specific = [ x for x in hg_j if ((x not in h1_j) and (x not in h2_j))]


    gtfStart, gtfEnd = self.readGTF()

    hap1_freq, hap1_oj,hap1_oj_count,hap1_novel = self.calculateFrequency(hap1_specific, h1_fn, gtfStart, gtfEnd)
    hap2_freq, hap2_oj,hap2_oj_count,hap2_novel = self.calculateFrequency(hap2_specific, h2_fn, gtfStart, gtfEnd)
    hap1hap2_freq1, hap1hap2_oj1,hap1hap2_oj_count1, hap1hap2_novel1 = self.calculateFrequency(hap1hap2_specific, h1_fn, gtfStart, gtfEnd)
    hap1hap2_freq2, hap1hap2_oj2,hap1hap2_oj_count2, hap1hap2_novel2 = self.calculateFrequency(hap1hap2_specific, h2_fn, gtfStart, gtfEnd)
    hg19_freq, hg19_oj,hg19_oj_count,hg19_novel = self.calculateFrequency(hg19_specific, hg_fn, gtfStart, gtfEnd)

    hap1hap2_freq = defaultdict(float)
    for j in hap1hap2_freq1:
        hap1hap2_freq[j] = (hap1hap2_freq1[j] + hap1hap2_freq2[j])/2

    hap1_out = os.path.join(self._outDir,"hap1.specific.bed")
    hap2_out = os.path.join(self._outDir,"hap2.specific.bed")
    hap1hap2_out = os.path.join(self._outDir,"hap1hap2.specific.bed")
    hg19_out = os.path.join(self._outDir,"hg19.specific.bed")

    self.printBed(hap1_specific,h1_j,hap1_freq,hap1_novel,hap1_out)
    self.printBed(hap2_specific,h2_j,hap2_freq,hap2_novel,hap2_out)
    self.printBed(hap1hap2_specific,h1_j,hap1hap2_freq,hap1hap2_novel1,hap1hap2_out)
    self.printBed(hg19_specific,hg_j,hg19_freq,hg19_novel,hg19_out)

    return


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
    outDir = open(".rPGAProject.yaml").readline().rstrip()
    vcf = open(".rPGAGenotype.yaml").readline().rstrip()
    gtf = open(".rPGAJunctions.yaml").readline().rstrip()
    ref = open(".rPGAGenome.yaml").readline().rstrip()
    seqs = open(".rPGASeqs.yaml")
    hap1Ref = os.path.join(outDir, "hap1.fa")
    hap2Ref = os.path.join(outDir, "hap2.fa")
    p = PersonalizeGenome(outDir, vcf, gtf, ref, hap1Ref, hap2Ref)
    if command == "personalize" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 1 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
      else :
        p.personalizeGenome()
        p.STAR_create_genome("HG19")
        p.STAR_create_genome("HAP1")
        p.STAR_create_genome("HAP2")
    elif command == "mapping" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 1 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
      else :
        p.STAR_perform_mapping("HG19", seqs)
        p.STAR_perform_mapping("HAP1", seqs)
        p.STAR_perform_mapping("HAP2", seqs)
    elif command == "discover" :
      if args[0].strip().lower() == "help" :
        print "Help"
      elif len(args) != 1 :
        sys.stderr.write("Genome file is not correct\n")
        sys.exit()
      else :
        p.haplotypeSpecificSam()
        p.distinctSJOut()
        p.haplotypeSpecificJunctions()
    else :
      sys.stderr.write("rPGA genomes -- unnknown command: " + command + "\n")
      sys.stderr.write(helpStr + "\n\n")
