BWA - GBS paired end files alignment to reference genome
===

This walkthrough shows how to align paired files to the reference (WGS) using BWA and checks how many reads have been mapped.
**Note** that this walkthrough assumes inputs as produces by the *GBS_raw2filtered* walkthrough.

Set up working space
In the desired directory, set up the working space required for this walkthrough:

```
mkdir align
```

This walkthrough starts in the 'paired' directory as created in the *GBS_raw2filtered* walkthrough.

A) Rename sequence identifiers of .fq files

During process.radtags, stacks renamed R1 and R2 sequence identifiers to _1 and _2. When running bwa mem it expects same file names or an error will appear:
[mem_sam_pe] paired reads have different names:
"5_1202_2403_86748_1", "5_1202_2403_86748_2"

Just cut of the _1 and _2 from the identifiers

```
sed 's/_1$//g' filenameR1 > newfilenameR1
sed 's/_2$//g' filenameR2 > newfilenameR2
for i in *.R1_f.fq; do basename=${i/.R1_f.fq}; sed 's/_1$//g' $i > ${basename}.R1.fq; done
for i in *.R2_f.fq; do basename=${i/.R2_f.fq}; sed 's/_2$//g' $i > ${basename}.R2.fq; done
```

**OR** (as I have encountered problems with above commandline:

```
ls -laht | awk '{print $NF}' | grep 'R1' | sort | uniq > list_allR1
ls -laht | awk '{print $NF}' | grep 'R2' | sort | uniq > list_allR2
while read i; do sed  's/_1$//g' $i > ${i}_s; done < list_allR1
while read i; do sed  's/_2$//g' $i > ${i}_s; done < list_allR2
```

Get all this to the align directory and rename the files

```
mv *fq_s ../align
rename .fq_s .fq *fq_s
```

---

B) Run bwa

This step includes conversion from .sam to .bam and sorting. Downstream we use steps and calling files without the "paired_" prefix, so this has to be removed from the files. Also, we want to lose the 'pl' from the names and call them 1-* instead of pl1-*...

```
for f in paired_pl*; do mv "$f" "${f#paired_pl}"; done
```


Make a list with the just the basenames (cut off last 6 characters, the .R1.fq and .R2.fq)
```
ls -laht | awk '{print $NF}' | grep '.fq$' | sed 's/\(.*\)....../\1/' | sort | uniq > list_sample
```

Load required modules

```
module load bwa
module load samtools
```

```
#bwa alignment script
#!/bin/bash
REF=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/ref/ragweed_10Feb2016_2ABsE_uppercase70.fasta
SAM_INDEX=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/ref/ragweed_10Feb2016_2ABsE_uppercase70.fasta.fai
OUT_PATH=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/output_all/bwa_genome_12567
#align to indexed database -t=threads -M=important for Picard -I = Illumina 1.3+ -q is quality trimming
bwa mem -t 2 -M $REF $1\.R1_tfq20.fq $1\.R2_tfq20.fq > $OUT_PATH/$1_tfq20.sam
#sam to bam conversion
samtools view -u -t  $SAM_INDEX -S $OUT_PATH/$1_tfq20.sam -o $OUT_PATH/$1_tfq20.bam
samtools sort $OUT_PATH/$1_tfq20.bam $OUT_PATH/$1_tfq20.sort
samtools index $OUT_PATH/$1_tfq20.sort.bam
```

And call the script with

```
while read i; do sh ~/scripts/SNP_calling/bwa_cra_lotte.sh $i; done < list_sample
```

Example output:

[M::main_mem] read 206186 sequences (20000042 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 30466, 47, 0)
[M::mem_pestat] skip orientation FF as there are not enough pairs
[M::mem_pestat] analyzing insert size distribution for orientation FR...
[M::mem_pestat] (25, 50, 75) percentile: (114, 166, 211)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 405)
[M::mem_pestat] mean and std.dev: (169.49, 71.38)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 502)
[M::mem_pestat] analyzing insert size distribution for orientation RF...
[M::mem_pestat] (25, 50, 75) percentile: (39, 39, 39)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (39, 39)
[M::mem_pestat] mean and std.dev: (39.00, 0.00)
[M::mem_pestat] low and high boundaries for proper pairs: (39, 39)
[M::mem_pestat] skip orientation RR as there are not enough pairs
[M::mem_pestat] skip orientation RF

This should be run in parallel and an example job is found in bwajob.job (align dir)


Check number of mapped and unmapped reads
the file produced will show the number of mapped reads, unmapped reads reads where
both read pairs are mapped and the bam file size (useful info on
http://left.subtree.org/2012/04/13/counting-the-number-of-reads-in-a-bam-file/)
the du command gives size of file.

````
echo -e 'sample\nmapped_reads\nunmapped_reads\nboth_reads_mapped\nbam_size' > mappedreads.txt;
while read i;  do echo "$i" >> mappedreads.txt; samtools view -c -F4 ${i}.sort.bam >> mappedreads.txt;  
samtools view -c -f 4 ${i}.sort.bam >>mappedreads.txt; samtools view -c -f1 -F12 ${i}.sort.bam >> mappedreads.txt; 
du ${i}.sort.bam | awk '{print  $1}' >> mappedreads.txt; done < list_sample

#and transform the rows to columns:
awk 'ORS=NR%5?" ":"\n"' mappedreads.txt > mappedreads.txt
```

