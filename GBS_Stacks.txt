===
STACKS

Can be used with or without a reference genome (this can be from a closely related species).
de-multiplexing

***STEP 1 - creating stacks that have a true polymorphism, allowing for difference due to mutation and sequencing error***
-identifying heterozygous loci within an individual (ustacks)
-reads within an individual are aligned with respect to one another (stacks) —> all reads in a stack are 100% identical
-stacks with depth below threshold (set by the user) will be set aside (secondary stacks)
-later on, secondary stacks are attempted to align to identified stacks
-stacks with > 2sd reads are exluded (lumberjack stacks) —> usually repetitive elements
-matching stacks of different individuals at all user-defined nucleotide distances, are aggregated together
-too many different stacks connected are discarded

Set by user
*depth threshold
*number of nucleotide differences with which to align stacks together (still in 1 barcoded individual)

Requirements of data:
+high average quality score of reads (>90%) or else reads will be discarded

***STEP 2 - inferring alleles***



Task X: Run "denovo.pl" This program creates stacks within each individual, then puts these stacks together to form catalog stacks (loci).

A) Create an environment to work in
in output_all
$ mkdir stacks
$ cd stacks
$ mkdir raw samples stacks paired assembled
$ ls
assembled  paired  raw	samples  stacks
Copy all the cleaned files into samples

B) Create a mysql database to import the results into. On the MCC, this requires using credentials to get access to mysql:
$ mysql --user=stacksuser --password=st4cks852 [command]

To create the database:

$ mysql --user=stacksuser --password=st4cks852 -e "CREATE DATABASE 125cra_radtags"

The database table definitions need to be send to the server to create all necessary components of the database, telling the database how to link tables.

$ mysql --user=stacksuser --password=st4cks852 temp_radtags < /opt/sw/stacks/1.13/share/stacks/sql/stacks.sql

To check that this worked:
```
$ mysql --user=stacksuser --password=st4cks852 
$ use temp_radtags;
$ show tables;
```

this should show a list of tables:
+------------------------+
| alleles                | 
| batches                | 
| catalog_alleles        | 
| catalog_annotations    | 
| catalog_genotypes      | 
| catalog_snps           | 
| catalog_tags           | 
| chr_index              | 
| fst                    | 
| genotype_corrections   | 
| markers                | 
| matches                | 
| pileup                 | 
| populations            | 
| ref_radome             | 
| samples                | 
| sequence               | 
| sequence_blast         | 
| snps                   | 
| sumstats               | 
| unique_tags            | 
+------------------------+
21 rows in set (0.00 sec)

To exit mysql:
```
$exit
```

C) Next, we can run denovo_map.pl

To print all the files in the required format for the denovo_map.pl type:
```
$ ls -l *.fq | awk '{print "-s", $9, "\\"}'
```
copy paste into the command and delete the last \

Remember that if you somehow cause this script to mess up, you will have to go back and re-create your mysql database and re-import the definitions again. This is because the database won't write over itself. Specify an output folder and create it before attempting to run denovo_map.pl

Parameters used:
-m:	Minimum number of raw reads required to make a stack, so the minimum depth
-M:	The maximum number of nucleotide mismatches expected between haplotypes at a locus within an individual. With high expected overall heterozygozity, this value can be set higher than default 3
-n:	The maximum number of nucleotide mismatches expected between any two haplotypes in the population. Follows same reasoning as -M
-T: 	The number of thread cores to run Stacks on. A higher number of cores used, means a lower time to complete analysis. Find out how many cores can be used

Multiple batches can be run in the one created database (as mentioned above). This can be useful when running the assembly with different parameters. Design the batches (make a list of all the different parameters to run) and for every batch only change the output file -o batch name -D and batch number -b.

```
denovo_map.pl -m 3 -M 3 -n 2 -T 2 -o /nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/output_all/stacks/stacks/batch1 -b 1 -B 125cra_radtags -D “125cra_batch1” \
-s AA-1-10.R1_f.fq \
-s AA-1-10.R2_f.fq \
-s AA-1-19.R1_f.fq \
-s AA-1-19.R2_f.fq \
-s AA-1-20.R1_f.fq \
-s AA-1-20.R2_f.fq \
-s AA-1-21.R1_f.fq \
-s AA-1-21.R2_f.fq \
-s AA-1-24.R1_f.fq \
-s AA-1-24.R2_f.fq \
-s AA-10-25.R1_f.fq \
-s AA-10-25.R2_f.fq \
-s AA-10-26.R1_f.fq \
-s AA-10-26.R2_f.fq \
-s AA-10-30.R1_f.fq \
-s AA-10-30.R2_f.fq \
-s AA-11-1.R1_f.fq \
-s AA-11-1.R2_f.fq \
-s AA-11-13.R1_f.fq \
-s AA-11-13.R2_f.fq \
-s AA-11-22.R1_f.fq \
-s AA-11-22.R2_f.fq \
-s AA-11-25.R1_f.fq \
-s AA-11-25.R2_f.fq \
-s AA-11-4.R1_f.fq \
-s AA-11-4.R2_f.fq \
-s AA-12-10.R1_f.fq \
-s AA-12-10.R2_f.fq \
-s AA-12-11.R1_f.fq \
-s AA-12-11.R2_f.fq \
-s AA-12-13.R1_f.fq \
-s AA-12-13.R2_f.fq \
-s AA-12-8.R1_f.fq \
-s AA-12-8.R2_f.fq \
-s AA-13-12.R1_f.fq \
-s AA-13-12.R2_f.fq \
-s AA-13-13.R1_f.fq \
-s AA-13-13.R2_f.fq \
-s AA-13-16.R1_f.fq \
-s AA-13-16.R2_f.fq \
-s AA-13-17.R1_f.fq \
-s AA-13-17.R2_f.fq \
-s AA-14-10.R1_f.fq \
-s AA-14-10.R2_f.fq \
-s AA-14-11.R1_f.fq \
-s AA-14-11.R2_f.fq \
etc
```

To browse the database in a firefox window, make sure you are logged in through the username@msgln6.its.monash.edu.au -Y and vlm001 -Y (-Y will give access to windows)
To open firefox:
```$ firefox &```
In the firefox window, go to: https://localhost/stacks/catalog.php?db=temp_radtags&id=1

--------------------------------
Task 3: Match the output from the denovo output to the Loblolly genome using the program BWA.

A) Convert stacks output to fasta format. First line starts with a ">" and has information about the run. Then there is a hard return to make a second line. The second line has the sequence. There are no quality scores associated with FASTA.

batch1: ```perl /data/GBS/catalogue2consensusfasta.pl clean_pineGBS/pineGBSstacks/batch_1.catalog.tags.tsv > batch1_pineGBSfull_consensus.fasta```

output: batch1_pineGBSfull_consensus.fasta

B) Input the consensus file that you created into BWA to align it to the loblolly genome.
 open byobu and change the working directory to clean_pineGBS/pineGBSstacks/

batch1: ```bwa mem -t 2 /data/Pine_genome/ptaeda.v1.01.fa.masked_cutN.noblank batch1_pineGBSfull_consensus.fasta > batch1_pineGBSfull_consensus.sam```