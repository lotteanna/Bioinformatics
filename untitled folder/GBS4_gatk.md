Use GATK to call variants
===

This task is divided into Task 3.1 using UnifiedGenotyper; and 3.2 using Haplotypecaller

Make list to add Read Groups (which will run in  gatk_cra_pseudo.sh, in step 3.1A) 
gatk is delicate and can't work without specified read groups. I made a list with
what I think are appropriate read-groups, with library being the illumina plate, platform
illumina, platform unit the barcodes and sample name the name of each individual.

---
Task 3.1 Calling SNPs and genotypes with UnifiedGenotyper
A) Re-align around the indels. 
Indels can cause alignment mismatches, especially near the end of an alignment and for
larger indels. This happens because in the alignment algorithm introducing an indel is
more 'expensive' than introducing a SNP. Therefore, the global alignment of all samples
needs to be considered. Reasons to re-align: a) known sites of indels; b) indels in the
original alignment (as can be seen in the CIGAR string); c) sites showing evidence of
indel (clustering of SNPs)
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

```
#SCRIPT SNP_calling/gatk_cra_pseudo.sh
#/bin/bash

REF=/nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa

RUN=$(echo $1 | awk -F"." '{print $1}')
LIB=$(echo $1 | awk -F"." '{print $3}')
#SAMPLE=$(echo $1 | awk -F"." '{print $2}')
ID=$(echo $1 | awk -F"." '{print $1"-"$2}')

echo $RUN $LIB $SAMPLE $ID

java -jar $PICARD/AddOrReplaceReadGroups.jar I=$ID.sort.bam O=$ID.rg.bam LB=$LIB PL=illumina PU=$RUN SM=$ID VALIDATION_STRINGENCY=SILENT  >> $ID.log 2>&1

#index
samtools sort $ID.rg.bam $ID.rg.sort
samtools index $ID.rg.sort.bam

#ID indels
java -Xmx15g -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -S LENIENT -R $REF -nt 10 -o $ID.bam.list -I $ID.rg.sort.bam > $ID.realign1.log

# realign around indels
java -Xmx10g -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar -T IndelRealigner -R  $REF -targetIntervals $ID.bam.list -I $ID.rg.sort.bam -o $ID.realign.bam > $ID.realign2.log
```

Implemented in gatkpseudo.job in bwa_genome_pseudo dir

```
#excluded from script
#while read f; do bash ~/scripts/SNP_calling/gatk_cra_pseudo.sh $f; done < listreadgroups.txt
```

Check if all files are present. If not (which happened to me), they can be checked by
using the following command (where list_sample2 contains ALL sample names:

```
ls -laht | awk '{print $9}' | grep 'rg.bam' | cut -f1 -d. > list_got
```

This will print list_got with all *rg.bam files in the working directory
```diff <(sort list_got) <(sort list_sample2)```

ND: USE GREP TO CHECK HOW MANY READS HAVE BEEN RE-ALIGNED

C) Base recallibration (look at bwa_red_gatk2.sh, started writing in gatk_cra2.sh)
 this step is only possible after a preliminary run that identifies SNPs and creates a vcf file

D) SNP and genotype calling with UnifiedGenotyper
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

With local realignment (implemented in hapl.job):

```
#SCRIPT hapl_cra.sh
#!/bin/bash

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-R /nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa \
-T HaplotypeCaller \
-I bam.list \
-o hapl.raw.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 20.0 \
-minPruning 3
```

```
sed -e 's/$/.rg.sort.bam/' list_sample2 > rgsortbam.list

```

```
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

---
Task 3.3. Calling haplotypes with Beagle
needs a vcf as input


                                                                                                                                  
---