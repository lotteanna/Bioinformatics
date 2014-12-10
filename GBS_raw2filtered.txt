Raw GBS paired-end reads to alignment input files
===

Important websites

http://creskolab.uoregon.edu/stacks/comp/process_radtags.php
https://genomicislands.wordpress.com/2013/11/20/how-to-perform-a-rad-seq-de-novo-assembly/
http://creskolab.uoregon.edu/stacks/pe_tut.php
http://www.simison.com/brian/Illumina_Rad_process_notes.html

---

**Task 1:** 
We need to rename all the reads to the sample that the barcode/lane matches to. 
Then we can put all the lanes together and know each sample name. The sequences are 
cleaned using the Illumina quality tags.

A) gzipped files (.gz) need to be unzipped. Do this by typing in the command line while in
directory of files that needed to be unzipped:

```
for f in *.gz; do gzip -d $f;done
```

B) All files are separated into number of the wells. These need to be concatenated into
1 single file for R1 (read 1) and R2 (backwards read without barcode, but is paired to R1
by number). To concatenate everything in the right order, make sure that the numbers in
the directory are correctly sorted. Example to type in comment line while in directory of 
files needed to be concatenated:
cat Ragweed6_NoIndex_L004_R1_*.fastq > Ragweed6_R1.fastq

C) process_radtags.sh (http://creskolab.uoregon.edu/stacks/comp/process_radtags.php)
Two ways to do this (I used option2). Either run every single file (001,002….) through 
process.radtags, or concatenate all the files as shown in B). 
Option1:

```
#script process.radtags.plate$.sh
#/bin/bash

# note by simon. if you have trouble decoding a bash script, you can run it in debug mode
# by boing /bin/bash -x
. /etc/profile
module load stacks

#In this example, we are dealing with paired-end reads. Normally Stacks should be able to
#deal with this using the -P flag, but I can’t get it working. Instead, we define the 
#pairs with $1, which requires an argument input in the command line, e.g. 001, for well 
#1. This mean that the script has to be run for every single well, and it will create a 
#subfolder for every single well (see below)

pair_1=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_NoIndex_L004_R1_$1.fastq.gz
pair_2=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_NoIndex_L004_R2_$1.fastq.gz

#process_radtags is only able to process barcodes of the same length. We have 6, 7, and 8 bp length barcodes in the same lane and need to feed this into the process using:
barcode_file=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_6.txt #path to barcode file
barcode_file2=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_7.txt
barcode_file3=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_8.txt
out_dir=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6_output/$1 #path to output

#make output directory, this is a subfolder for every single well
mkdir $out_dir

#run process_radtags, in this example we have paired-end, Illumina HiSeq data and use the 
#-P flag to show this to the program. Also, data are gzipped, and the -i flag is specified:
#-i gzfastq

process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file -e pstI -c -q -r -i gzfastq
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file2 -e pstI -c -q -r -i gzfastq
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file3 -e pstI -c -q -r -i gzfastq
```

Option 2:

```
#script process.radtags.plate*all.sh

#/bin/bash
#see option 1 for better descriptions
. /etc/profile
module load stacks
pair_1=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_R1.fastq
pair_2=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_R2.fastq
barcode_file=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_6.txt #path to barcode file
barcode_file2=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_7.txt
barcode_file3=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_8.txt
out_dir=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6_output/all #path to output
mkdir $out_dir
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file -e pstI -c -q -r
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file2 -e pstI -c -q -r
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file3 -e pstI -c -q -r
```

to run the script: type ```sh process.radtags.plate[whichever you wanna run][all].sh
[well# whichever you wanna run, in case of option1, e.g. 004]```

These process_radtag output files are in ~/ragweed/GBS/raw_common/output_all/plate*_pr 
(pr for process radtags

D) Use rename.sh to rename the barcodes. Input the file that contains both the barcode 
and the sample name and this shell script will change the names.

```
#script rename.sh
#!bin/bash

#awk -F specifies the field separator, in this case .
barcode=$(echo $1 | awk -F"." '{print $2}')
sample=$(echo $1 | awk -F"." '{print $1} ')

echo $barcode
echo $sample

mv sample_$barcode.1.fq $sample.R1.fq
mv sample_$barcode.2.fq $sample.R2.fq
```

Then, type in command line (in directory where files have to be renamed):

```
for f in `cat /nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes_plate1.txt`; 
do sh ~/scripts/rename.sh $f; done
```

NB: this is to rename file names of plate1

Good practise to copy files into another folder so original process.radtags output files 
are saved. Renamed files (including .rem.fq files) are in 
~/ragweed/GBS/raw_common/output_all/plate*_rn

E) Next, we move all fq files into one folder. Create a new folder and move the 
renamed files into it. E.g. ```cp -r *R1.fq ../allind```

F) We need to remove reads that contain the adapter sequence or part of the adapter sequence.
Stacks is able to do this, but does a shit job. Kay made a perl script which is better, 
adapter_removal.pl

```
#!bin/perl
use warnings;
use strict;
#get file
my $fq_1=$ARGV[0];

open FQ1, $fq_1;

open FQ1_clean, ">$fq_1\_clean";
my $i=0;
my $header1 = ();
my $read1 = ();
my $qual1 = ();
my $badreadno1 =0;
while (<FQ1>) {
    my $line1 = $_;
    chomp $line1;
        if($i==0){
                if($line1=~m/^\@/){
                        $header1=$line1;
                        $i=1;
                }elsif($line1=~m/^\+$/){
                        $i=2;
                }
        }elsif($i==1){
                $read1=$line1;
                $i=0;
        }elsif($i==2){
                $qual1=$line1;
                $i=0;
                #remove forward adapter contamination whole or parts
                if($read1=~m/CTGCAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG/ || $read1=~m/CTGCAAGATCGGAAGAGCG/ || $read1=~m/CGGTTCAGCAGGAATGCCGAG/ || $read1=~m/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT/ || $read1=~m/AGATCGGAAGAGCGTCGTGT/ || $read1=~m/AGGGAAAGAGTGT/){
                        $badreadno1 = $badreadno1 +1;
                }else{
                        print FQ1_clean "$header1\n$read1\n+\n$qual1\n";
                }
        }

}
close FQ1_clean;
close FQ1;

print "The number of reads with adapter contamination is $badreadno1\n";
```
                                                                                
In allind type in command line:

```
for i in *.fq; do perl ~/scripts/adapter_removal.pl $i ; done
```

When working, this will give output lines on the screen:

>The number of reads with adapter contamination is 14329
>The number of reads with adapter contamination is 10637
>The number of reads with adapter contamination is 163986
>The number of reads with adapter contamination is 129353
>The number of reads with adapter contamination is 35568
>The number of reads with adapter contamination is 28743

This will additionally change the *.fq to *.fq_clean, so you know that the file has been
looked at. Check and compare file size to make sure files were reduced in size.

Copy the cleaned, renamed files into another new folder (clean):

```
cp -R *fq_clean ../clean
```

The output files have the ending '.fq_clean' Other programs do not like this. Change the
endings by using the command:

```
rename .fq_clean .fq *fq_clean
```

E) Trimming reads.  NOTE: This step has to be skipped when using STACKS, as it is not
able to handle reads of different lengths. bwa is though.
Sickle uses sliding windows along with quality and length thresholds
to determine when quality is sufficiently low to trim both ends of reads.
the following option in Sickle takes 2 input files, produces 2 output files and a
orphan file:

```
sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
-o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
-s trimmed_singles_file.fastq
```

pe stands for paired-end reads, -t is quality encoding (phred+33 is sanger)

```
for i in *R1.fq; do sickle pe -t sanger -f $i -r ${i/R1.fq/}\R2.fq \
-o ../trimmer/${i/R1.fq/}\R1_tr.fq -p ../trimmer/${i/R1.fq/}\R2_tr.fq \
-s ../trimmer/${i/R1.fq/}\S_tr.fq -q 20 -l 50 \
>> ../trimmer/${i/.R1.fq/}\_tr.log 2>&1 $i; done
```

F) It is really important that all the reads are of good quality because STACKS aligns
them without caring about the quality score. This step might be skipped when trimmer is
already used, but I prefer to use both, as the trimmer only trims off from the ends
until the quality is good. But maybe the average read qual is still crap
load fastx
. /etc/profile
module load fastax
run quality filter on .fq.clean and creates output files _q10.fq and .fq.log (log file
for filter). Or _tfq20.fq for both trimmed and filtered reads (in dir trimfilt)

```
for i in *.fq; do  fastq_quality_filter -q 10 -p 90 -Q 33 -v -i $i \
-o ../filter/q10/${i/.fq/}\_q10.fq >> $i_q10.log 2>&1 $i; done

for i in *.fq; do  fastq_quality_filter -q 20 -p 90 -Q 33 -v -i $i \
-o ../filter/q20/${i/.fq/}\_q20.fq >> $i_q20.log 2>&1 $i; done
```

To see the amount of reads that are filtered out, I can used the produced log files
first make a list with the names

```
ls -laht | awk '{print $NF}' | grep ‘_q10.fq$’ | cut -d. -f1  | sort | uniq > list_sample
awk '{print $0".R1"}' list_sample > list_sampleR1R2
awk '{print $0".R2"}' list_sample >> list_sampleR1R2
```

Then make a new text file containing this information

```
echo -e "Sample\nInput\nOutput\nDiscarded\n" > qualfiterlog.txt; 
while read i;  
do echo "$i" >> qualfiterlog.txt; 
cat ${i}.fq.clean.log | awk '$1 ~ /^Input:/ { print $2}' >> qualfiterlog.txt;  
cat ${i}.fq.clean.log | awk '$1 ~ /^Output:/ { print $2}' >> qualfiterlog.txt; 
cat ${i}.fq.clean.log | awk '$1 ~ /^discarded/ { print $2}' >> qualfiterlog.txt; 
done < list_sampleR1R2
```

Remove whitespaces and empty lines respectively:

```
sed -i -e 's/^[ \t]*//' -e 's/[ \t]*$//' qualfiterlog.txt
sed -i '/^$/d' qualfiterlog.txt
```

Put all this in 4 columns

```
awk 'ORS=NR%4?" ":"\n"' qualfiterlog.txt >> qualfiterlog2.txt
```

Copy to local computer (not ssh into mcc)

```
scp -pr lotteanv@msgln4.its.monash.edu:~/ragweed/GBS/raw_common/output_all/clean/qualfiterlog2.txt Documents
```

G) Re-pair reads. During filtering some reads have been filtered out, which means that
gaps/differences exist between R1 and R2. This can be fixed by re-pairing the R1 and R2
files. NOTE however that below script does not do anything with "orphans", non-paired reads

Copy all _f.fq files to pairing:

```cp -R *_f.fq ../pairing```

Saved in ~scripts/fix_fqpair_all.sh
Let's go to your directory where the files live
cd ~/ragweed/GBS/raw_common/output_all/pairing

Throw all basenames of files to be processed into list. here we look at all files ending
with .fq - change your grep if you want to be more specific. make sure you dont capture
out output files into list_sample. also this assumes the base file name is before the
first dot.

```
ls -laht | awk '{print $NF}' | grep '_q10.fq$' | cut -d. -f1 | sort | uniq > list_sample
```

Do shit to the list.

```
while read i
do
perl ~/scripts/FQ_pair_no.pl ${i}.R1_q10.fq ${i}.R2_q10.fq >> $i_q10.log 2>&1
done < list_sample
```

```
mv paired* ../paired
```







