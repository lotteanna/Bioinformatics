Use GATK to call variants
===

This task is divided into Task 3.1 using UnifiedGenotyper; and 3.2 using Haplotypecaller

Make list to add Read Groups (which will run in  gatk_cra_pseudo.sh, in step 3.1A) 
gatk is delicate and can't work without specified read groups. I made a list with
what I think are appropriate read-groups, with library being the illumina plate, platform
illumina, platform unit the barcodes and sample name the name of each individual.

---
Calling SNPs and genotypes with UnifiedGenotyper

**A) Re-align around the indels **

Indels can cause alignment mismatches, especially near the end of an alignment and for
larger indels. This happens because in the alignment algorithm introducing an indel is
more 'expensive' than introducing a SNP. Therefore, the global alignment of all samples
needs to be considered. Reasons to re-align: 

a) known sites of indels; 

b) indels in the
original alignment (as can be seen in the CIGAR string); 

c) sites showing evidence ofindel (clustering of SNPs)

It is good to include sites of known indels & sites of known variation in the realignment.
RealignTargetCreator creates targets of interest for the realignment and outputs an
.intervals file. This tool can use a vcf file of known indels. Not necessary if not
available, but will speed up process.
The IndelRealigner is a heavier computational process which is doing the actual realignment
It changes the CIGAR string of re-aligned reads and maintains the original CIGAR string
with an OC tag. It is thus easy to assess afterwards how many reads have been re-aligned
by using grep etc.
A full Smith-Waterman realignment is recommended when: a) older data; b) shorter reads
(~36bp); c) no known indels; d) if you want to reduce false-positives
Indel realignment is very important before proceeding with the Base Quality Score
Recalibration

---

First step is to make sure everything is copied in the correct directories and rename the sample files to the core sample name (as in listreadgroup.txt, below):

```
mkdir ../realign
cp * ../realign
cd ../realign
rename _tfq20. . *_tfq20.*
```

Secondly, you need to make a text file called listreadgroups.txt. This file will contain information about the plate, the sample name and some identifier. I search hard and long what this actually means (as it will be used in the gatk pipeline), but I couldn't get any further then specifying the lane number, sampleID and barcode used. So that is exactly listreadgroups,txt, separated by a ".":

1.010908-1-27.CTCGCA
1.010908-1-12.CCAGCT
1.160808-1-12.TAGGTCA
1.160808-1-3.TAATTCG
1.160808-1-1.TTCAGAG
1.160808-2-9.CCAGGTA
1.160808-2-8.CAAGCTG
1.160808-2-7.AATAAGCG

etc



Now, the below script will use this information to add something to the .rg.bam output files using the Picard function AddOrReplaceReadGroups. Also, it will index the produced files, create targets for the indel realignment and finally, realign around the indels.

The gatk pipeline is a little picky, and doesn't work with all versions of everything. I've tried a few things, and this seems to work (you need to start a new terminal window to be able to implement these):

```
module load java/jdk1.7.0_21
module load gatk/3.1.1
module load samtools
```

And the script

```
#SCRIPT SNP_calling/gatk_cra_full.sh
#/bin/bash

REF=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/ref/ragweed_10Feb2016_2ABsE_uppercase70.fasta.fa
RUN=$(echo $1 | awk -F"." '{print $1}')
LIB=$(echo $1 | awk -F"." '{print $3}')

#SAMPLE=$(echo $1 | awk -F"." '{print $2}')
ID=$(echo $1 | awk -F"." '{print $1"-"$2}')
echo $RUN $LIB $SAMPLE $ID
java -jar /opt/sw/picard/1.1.08/picard-tools-1.108/AddOrReplaceReadGroups.jar I=$ID.sort.bam O=$ID.rg.bam LB=$LIB PL=illumina PU=$RUN SM=$ID VALIDATION_STRINGENCY=SILENT  >> $ID.log 2>&1

#index
samtools sort $ID.rg.bam $ID.rg.sort
samtools index $ID.rg.sort.bam

#ID indels
java -Xmx15g -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -S LENIENT -R $REF -nt 10 -o $ID.bam.list -I $ID.rg.sort.bam > $ID.realign1.log

# realign around indels
java -Xmx10g -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar -T IndelRealigner -R  $REF -targetIntervals $ID.bam.list -I $ID.rg.sort.bam -o $ID.realign.bam > $ID.realign2.log
```

Implemented in gatkpseudo.job in bwa_genome dir


Check if all files are present. If not (which happened to me), they can be checked by
using the following command (where list_sample2 contains ALL sample names:

```
ls -laht | awk '{print $9}' | grep 'rg.bam' | cut -f1 -d. > list_got
```

This will print list_got with all *rg.bam files in the working directory
```diff <(sort list_got) <(sort list_sample2)```

Move all files to gatk_12567

ND: USE GREP TO CHECK HOW MANY READS HAVE BEEN RE-ALIGNED

C) Base recallibration (look at bwa_red_gatk2.sh, started writing in gatk_cra2.sh)
 this step is only possible after a preliminary run that identifies SNPs and creates a vcf file

D) SNP and genotype calling with UnifiedGenotyper
UG is the main tool to get genetic variance. Most of the time, real mutations are hidden in the noise; need analysis tools to determine if variation is genetic variance of random machine noise. UG is Bayesian and calls SNPs and indels seperatly by considering each variation independently and assigns a genotype to a sample when the variant is called.

Make a list of all the realigned bam files to be the input file for the unified genotyper

```
ls -laht | awk '{print $9}' | grep 'realign.bam' > bam.list
```

```
#script gatk_cra3.sh
#!/bin/bash

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-R /nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa \
-T UnifiedGenotyper \
-I bam.list \
-o snps.raw.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 20.0 \
```

```sh gatk_cra3.sh```

when logging a job for 288 samples, this took 2.5 days
---
Task 3.2. Calling haplotypes with HaplotypeCaller
This is a fairly new tool and works at the moment only with diploid organisms
Approximately 5 times slower than UnifiedGenotyper

A) Haplotype calling
HaplotypeCaller is approximately 5 times slower than UnifiedGenotyper, depending on
parameters set by user. 
This tool is able to detect large variance and looks at regions of variation instead of 
independent loci. It does this through local denovo assembly. It uses all the reads
localised to a region, so it does not assume the exact alignment made by the alignment
software. For this reason, indel re-alignment is not necessary (but it should not matter
if this step is preformed before using this tool). HC does not use prior known data,
like known indels etc.
This tool works better than UG because it takes multiple loci into account, e.g UG
calls several loci where it's only 1 deletion (same reason for ID realignment). HC also
calls variants on regions that are too big for ID realigner.
HC brings in ALL mate pairs, including un-mapped mates.

Make list of input bam files (I WANT TO CHECK IF RESULTS ARE DIFFERENT WITH AND WITHOUT
LOCAL INDEL REALIGNMENT)


```
ls -laht | awk '{print $9}' | grep 'rg.sort.bam'  > rgsortbambai.list
cat rgsortbambai.list | grep -v '.bai' > rgsortbam.list
```

With local realignment (implemented in hapl.job):

```
#SCRIPT hapl_cra.sh
#!/bin/bash

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-R /nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa \
-T HaplotypeCaller \
-I rgsortbam.list \
-o hapl.raw.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 20.0 \
-minPruning 3

```

Pruning: the amount of minimum reads that have to show the variation to be included.

Interesting read on:
WARN  22:20:14,941 ExactAFCalc - this tool is currently set to genotype at most 6 alternate alleles in a given context, but the context at Scaffold_586:494835 has 7 alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument

                                                                                                                                  
---