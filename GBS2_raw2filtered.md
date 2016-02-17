Raw GBS paired-end reads to alignment input files
===

Load required modules
```
module load fastqc
module load sickle
```

This walkthrough shows how rename and filter raw paired end Illumina. Firstly, all the reads need to be renamed to the sample that the barcode/lane matches to – *de-multiplexing. Then, we will remove adaptors, trim, filter and re-pair the reads. 

Note that this walkthrough is followed by the *GBS_bwa* walkthrough, and produces files that are referred to downstream

——
Important websites

http://creskolab.uoregon.edu/stacks/comp/process_radtags.php
https://genomicislands.wordpress.com/2013/11/20/how-to-perform-a-rad-seq-de-novo-assembly/
http://creskolab.uoregon.edu/stacks/pe_tut.php
http://www.simison.com/brian/Illumina_Rad_process_notes.html
http://www.ark-genomics.org/events-online-training-eu-training-course/adapter-and-quality-trimming-illumina-data

---

Set up the working space. In a directory where you want all the output files, make the required output directories
```
mkdir plate[plateID]_pr plate[plateID]_rn clean trimmed filtered pairing paired
```

**A) Unzip files**
Gzipped files (.gz) need to be unzipped. Do this by typing in the command line while in
directory of files that needed to be unzipped

```
for f in *.gz; do gzip -d $f;done
```

**B) Concatenate files per lane**
All files are separated into number of the wells. These need to be concatenated into
1 single file for R1 (read 1) and R2 (backwards read without barcode, but is paired to R1
by number). To concatenate everything in the right order, make sure that the numbers in
the directory are correctly sorted. 

```
cat Ragweed6_NoIndex_L004_R1_*.fastq > Ragweed6_R1.fastq
```

**C) Explore quality of reads**
This is a perfect time to explore the quality of the raw sequences. All that you need to do is

```
fastqc *.fastq
```

This will produce a html file with a lot of information that can be used for later filtering steps (and the number of reads in the raw data!) You need to download the file or use a browser in your terminal
This step can be repeated after every filtering step to see how the process is going. 

**D) De-multiplex reads**
The concatenated R1 and R2 files are a soup of all the barcoded individuals and their reads. We need to ‘de-multiplex’ the reads to get all the individuals and their corresponding reads, based on their barcodes. This can be done with: process_radtags.sh (http://creskolab.uoregon.edu/stacks/comp/process_radtags.php)

I made a script for every single plate, as all the infiles were in different directories. It can be done in one go though. Below script is for plate 6.

```
#script process.radtags.plate*all.sh
#/bin/bash
. /etc/profile
module load stacks
pair_1=/<dir>/Ragweed6_R1.fastq
pair_2=/<dir>/Ragweed6_R2.fastq
#process_radtags is only able to process barcodes of the same length. We have 6, 7, and 8 bp length barcodes in the same lane and need to feed this into the process using:
barcode_file=/<dir>/barcodes97to192_6.txt #path to barcode file
barcode_file2=/<dir>/barcodes97to192_7.txt
barcode_file3=/<dir>/barcodes97to192_8.txt
out_dir=/<dir>/plate6_pr #path to output
mkdir $out_dir
#run process_radtags, in this example we have paired-end, Illumina HiSeq data and use the 
#-P flag to show this to the program. Also, data are gzipped, and the -i flag is specified:
#-i gzfastq
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file -e pstI -c -q -r
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file2 -e pstI -c -q -r
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file3 -e pstI -c -q -r
```

Run the script

```
sh process.radtags.plate[plate#]all.sh
```

**E) Rename the barcodes**
The barcodes don’t hold much information. We want the sample names. So let’s rename all samples. Input the file that contains both the barcode and the sample name and this shell script will change the names.

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

Go to plate[plate#]_pr and call the script

```
for f in `cat /nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes_plate1.txt`; 
do sh ~/scripts/rename.sh $f; done
```

NB: this is to rename file names of plate1

Good practise to copy files into another folder so original process.radtags output files 
are saved. 

```
mv sample* ../plate[plate#]_rn
```

**F) Copy files to filter dirs**
Next, we copy all fq files into a new dir (so we won’t mess up our precious renamed files). 

```
cp -r *R1.fq ../allind
cp -r *R2.fq ../allind
```

**G) Remove adaptor sequences**
We need to remove reads that contain the adapter sequence or part of the adapter sequence. Multiple packages exist to do this, including Scythe and Stacks. The latter however does shit job. A better custom script is below, which will actually specify the sequences of the adaptor, and takes them out.

```
#adapter_removal.pl
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
                                                                                
Call the script

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

This will additionally change the *.fq to *.fq_clean, so you know that the file has been looked at. Check and compare file size to make sure files were reduced in size.

Copy the cleaned, renamed files into another new folder (clean) and rename them.

```
cp -R *fq_clean ../clean
rename .fq_clean .fq *fq_clean
```

H) Trimming reads  

```mkdir trimmer```

*NOTE: This step has to be skipped when using STACKS, as it is not able to handle reads of different lengths. bwa is though.*

The ends of reads are always of lower quality, and that is why we would like to take of a little bit of each end to make the overall read score better. You can also explore this in your fastqc plots, and can use these plots to determine how much you should trim of each end.

One package to use is sickle (https://github.com/najoshi/sickle).
Sickle uses sliding windows along with quality and length thresholds to determine when quality is sufficiently low to trim both ends of reads. 

There are several options in Sickle, including ones that deal with a single 'interleaved' paired end input files, or 2 seperate input files for both paired-end reads.
The following option in Sickle takes 2 paired-end input files, produces 2 output files and an "singles" file. The "singles" file containes reads that passed the sickle trimmer in either the forward or reverse direction, but not both.

```
sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
-o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
-s trimmed_singles_file.fastq
```

pe stands for paired-end reads, -t is quality encoding (phred+33 is sanger)

And now a real example (it's possible below code won't work directly due to the \. Just delete these and write as one sentence)

```
for i in *R1.fq; do sickle pe -t sanger -f $i -r ${i/R1.fq/}\R2.fq \
-o ../trimmer/${i/R1.fq/}\R1_tr.fq -p ../trimmer/${i/R1.fq/}\R2_tr.fq \
-s ../trimmer/${i/R1.fq/}\S_tr.fq -q 20 -l 50 \
>> ../trimmer/${i/.R1.fq/}\_tr.log 2>&1 $i; done
```

Now, check out the log file. I had an interesting warning:
```Warning: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.```


F) Quality filtering
It is really important that all the reads are of good quality because STACKS aligns them without caring about the quality score. This step might be skipped when trimmer is already used, but I prefer to use both, as the trimmer only trims off from the ends until the quality is good. But maybe the average read quality is still crap.

Load fastx

```
. /etc/profile
module load fastax
```

run quality filter on .fq.clean and creates output files _q10.fq and .fq.log (log file for filter). Or _tfq20.fq for both trimmed and filtered reads 

```
for i in *.fq; do  fastq_quality_filter -q 10 -p 90 -Q 33 -v -i $i \
-o ../filter/q10/${i/.fq/}\_q10.fq >> $i_q10.log 2>&1 $i; done
for i in *.fq; do  fastq_quality_filter -q 20 -p 90 -Q 33 -v -i $i \
-o ../filter/q20/${i/.fq/}\_q20.fq >> $i_q20.log 2>&1 $i; done
for i in *.fq; do  fastq_quality_filter -q 20 -p 90 -Q 33 -v -i $i -o \
${i/_tr.fq/}\_tfq20.fq >> ${i/_tr.fq/}\_tfq20.log 2>&1 $i; done
```

Notice that above code will write all log files to the trimmer folder (as the output wasn't set to another folder. This can be fixed by either moving the log files OR (code not tried):

```
for i in *.fq; do  fastq_quality_filter -q 20 -p 90 -Q 33 -v -i $i -o \
../filter/${i/_tr.fq/}\_tfq20.fq >> ../filter/${i/_tr.fq/}\_tfq20.log 2>&1 $i; done
```

To see the amount of reads that are filtered out, I can used the produced log files first make a list with the names.

```
ls -laht | awk '{print $NF}' | grep ‘_q10.fq$’ | cut -d. -f1  | sort | uniq > list_sample
awk '{print $0".R1"}' list_sample > list_sampleR1R2
awk '{print $0".R2"}' list_sample >> list_sampleR1R2
```

Then make a new text file containing this information.

```
echo -e "Sample\nInput\nOutput\nDiscarded\n" > qualfiterlog.txt; 
while read i;  
do echo "$i" >> qualfiterlog.txt; 
cat ${i}.fq.clean.log | awk '$1 ~ /^Input:/ { print $2}' >> qualfiterlog.txt;  
cat ${i}.fq.clean.log | awk '$1 ~ /^Output:/ { print $2}' >> qualfiterlog.txt; 
cat ${i}.fq.clean.log | awk '$1 ~ /^discarded/ { print $2}' >> qualfiterlog.txt; 
done < list_sampleR1R2
```

Remove white spaces and empty lines.

```
sed -i -e 's/^[ \t]*//' -e 's/[ \t]*$//' qualfiterlog.txt
sed -i '/^$/d' qualfiterlog.txt
```

Put all this in 4 columns.

```
awk 'ORS=NR%4?" ":"\n"' qualfiterlog.txt >> qualfiterlog2.txt
```

Copy to local computer (not ssh into mcc)

```
scp -pr lotteanv@msgln4.its.monash.edu:~/ragweed/GBS/raw_common/output_all/clean/qualfiterlog2.txt Documents
```

I) Re-pair reads 
During filtering some reads have been filtered out, which means that gaps/differences exist between R1 and R2. This can be fixed by re-pairing the R1 and R2 files. NOTE however that below script does not do anything with "orphans", non-paired reads.

Copy all non-singleton files to pairing:

```
cp -R *.R* ../pairing
```

Let's go to your directory where the files live

```
cd ../pairing
```

Throw all basenames of files to be processed into list. Here we look at all files ending with .fq - change your grep if you want to be more specific. Make sure you don’t capture
out output files into list_sample. also this assumes the base file name is before the first dot.

```
ls -laht | awk '{print $NF}' | grep '.fq$' | cut -d. -f1 | sort | uniq > list_sample
```

Do shit to the list.

```
while read i; do perl ~/scripts/FQ_pair_no.pl ${i}.R1_tfq20.fq ${i}.R2_tfq20.fq >> $i_tfq20.log 2>&1; done < list_sample
```

```
mv paired* ../paired
```
—

Continue with the GBS_bwa walkthrough







