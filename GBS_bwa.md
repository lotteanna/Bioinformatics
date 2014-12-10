BWA - GBS paired end files alignment to fragmented reference genome

===

Align paired files to crappy reference (WGS) using BWA (NB: I chopped everything up in
one-liners but these can be put together in a script)
```
. /etc/profile
module load bwa
```

A) Rename sequence identifier of fq files. 

During process.radtags, stacks renamed R1 and R2 sequence identifiers to _1 and _2. When running bwa mem it expects same file names or an error will appear, e.g. [mem_sam_pe] paired reads have different names: 
"5_1202_2403_86748_1", "5_1202_2403_86748_2" 
```sed 's/_1$//g' filenameR1 > newfilenameR1```
```sed 's/_2$//g' filenameR2 > newfilenameR2```

```for i in *.R1_f.fq; do basename=${i/.R1_f.fq}; sed 's/_1$//g' $i > ${basename}.R1.fq; done```
```for i in *.R2_f.fq; do basename=${i/.R2_f.fq}; sed 's/_2$//g' $i > ${basename}.R2.fq; done```

OR (as I have encountered problems with above commandline:

```
ls -laht | awk '{print $NF}' | grep 'R1’ | sort | uniq > list_allR1
ls -laht | awk '{print $NF}' | grep 'R2' | sort | uniq > list_allR2
while read i; do sed  ’s/_1$//g' $i > ${i}_s; done < list_allR1
while read i; do sed  's/_2$//g' $i > ${i}_s; done < list_allR2
```

Copy *fq_s to align folder

```
rename .fq_s .fq *fq_s
```

B) (optional) In case of fragmented reference, it is better to make a pseudo-scaffold. This is 
a concatenated file with all the contigs togethers, separated with 30 A's.
perl ~/scripts/scaffolds.pl soaprunk61.contig
output files are soaprunk61.contig.pseudo  & scaffold_order.soaprunk61.contig
REMEMBER to translate the SNP table back later to the original genome locations
 with scaffold2contig.v2.pl

C) GATK is not able to handle N's, count and replace these with A's
Count number of N’s in the file as GATK can’t handle N’s
note: whenever counting lines make sure to account for header lines.
Use grep -v '>'  first bascially because you don’t want to count characters in the
header and you should never ever assume the header won't contain some chars.
grep -v is "inverse grep": that is "look for lines without this"
But without accounting for this:

**Option1** (—> gives string of A’s and counts the characters)
tr -cd N < soaprunk61.contig.pseudo | wc -c 
 
*OR*

**Option2**(—> slower option, gives char in lines and then counts the number of lines.
 wc -c won’t work here as there will be twice as many characters, A and newline)
```fgrep -o N soaprunk61.contig.pseudo | wc -l```

In case of N's replace with A's with sed, not done here as wasn't necessary

can also be used to see how many A's were added in the scaffolds.pl script
by counting the difference between soaprunk61.contig and soaprunk61.contig.pseudo
```echo '712102436 - 421159706' | bc```

D) Index .contig (or in case of B, contig.pseudo) file for samtools, bwa and gatk

Check maximum line lengths in file with:
```awk '{ if (x < length()) x = length() } END { print x }' <referencegenome>```
Check for anything shorter that this maximum length (here 70)
```awk 'length($0) < 70' <referencegenome>```
In case not all lines are same length:

Delete white lines
```sed '/^$/d' <filename> > <newfilename>```

-a specifies the indexing algorithm bwa uses. There is another option (IS), for smaller
genomes (<2GB), check bwa man page
```bwa index -a bwtsw <referencegenome>```

```samtools faidx <referencegenome>```

Note by Kay (SNP calling blog): if you want to use Picard downstream use the -M option in bwa. This is NOT
implemented in below script! From bwa man: Mark shorter split hits as secondary


E) Run bwa (in align folder) 
this includes conversion from .sam to .bam and sorting
The listreadgroups list used later is without the "paired_" prefix, so this has to be
removed from the files. Also, this list has 1.* instead of pl1-* etc
```for f in paired_pl*; do mv "$f" "${f#paired_pl}"; done```
```rename pl1- 1- pl1-*```

etc

Make a list with the just the basenames (cut off last 6 characters)
```ls -laht | awk '{print $NF}' | grep '.fq$' | sed 's/\(.*\)....../\1/' | sort | uniq > list_sample```

```
#excluded from script
#while read i; do sh ~/scripts/SNP_calling/bwa_cra_lotte.sh $i; done < list_sample
#The following is implemented in bwapseudo.job in align directory:
#while read i; do sh ~/scripts/SNP_calling/bwa_cra_pseudo.sh $i; done < list_sample
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

To go through conversion of sam to bam as one-liners (this is however implemented in
bwa_cra_pseudo.sh script, not in bwa_cra_lotte.sh script):

```
#excluded from script
#for i in *.sam; do basename=${i/.sam};  samtools view -u -t  $SAM_INDEX  -S ${basename}.sam -o ${basename}.bam; done
#for i in *.bam; do basename=${i/.bam}; samtools sort ${basename}.bam ${basename}.sort; done
#for i in *sort.bam; do basename=${i/.sort.bam}; samtools index ${basename}.sort.bam; done
```

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

