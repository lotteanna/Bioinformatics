===
Task 4: Prepare the SNP-tables for the different programs

**STRUCTURE**
http://pritchardlab.stanford.edu/structure.html

A)

STRUCTURE can't deal with linked loci, so need to select at random one SNP per contig

```
perl /nfs/home/hpcsci/lotteanv/scripts/random_one_per_locus_combined.pl snptableUG.tab.table.contig
```

Check if number of selected contigs is the same as number of unique contigs
```
awk 'NR!=1{print $1}' snptableUG.tab.table.contig | sort | uniq -c | wc -l
```

5244
```
awk 'NR!=1{print $1}' snptableUG.tab.table.contig.random | wc -l
```

B)
Convert the SNPtable to a STRUCTURE readable format

 NOTE: the order of the digital numbers is mixed up (eg AC can now be read in by Structure as CA. This means linkage mapping is impossible with this data (see Structure documentation). This is caused by the coding used initially, which is by vcf2vertical_dep_GATK_UG.pl

First, make a population document (popind) for the conversion script, which is a tab separated file. The following seems to work (don't enter the headers):
pop        sampleID
010908-1	2-010908-1-17
010908-1	2-010908-1-9
010908-1	1-010908-1-27
010908-1	1-010908-1-12
010908-1	6-010908-1-7
010908-1	6-010908-1-22
160808-1	2-160808-1-11
160808-1	1-160808-1-12
etc

The original file I have to translate is called snptableUG.tab.192.table.contig.random. This file is structured as follows:

Columns:
1: Chrom#
2: Locus Location
3 - 386: SamplesIDs

This table doesn't neatly have the genotypes spelled out as AA, AT, GC etc, but uses the combined coding for genotypes: A, N, W, etc.

First of all, as I understand from the STRUCTURE manual and your 'structure_input_file.R' script, I have to divide the data into 2 rows per individual, each row for one of the 2 alleles (as ragweed is diploid). So my plan is to:
- Digitalise and separated genotypes for each allele ('A' = 1,1; 'N' = -9,-9; 'W' = 1,2 etc)
- Concatenate Chrom and Locus Location
- Add in population information for each individual (through popind file)
- Transpose columns and put them in the right order for STRUCTURE

Run the conversion script:
```
#########################################################################
# SNPtable2structure_combined_lotte                                     #
#                                                                       #
# when including location information run as:                           #
# perl bin/SNPtable2structure.pl data/rand.snptable \                   #
# results/rand.structure data/names2pop                                 #
# when not including location information run as:                       #
# perl bin/SNPtable2structure.pl data/rand.snptable \                   #
# results/rand.structure                                                #
#                                                                       #
# This script converts custom SNP tables to STRUCTURE format            #
# Current version allows for 2 SNPlocation information columns          #
# last version: transcript contig, genomic contig and genomic position  #
# current version: genomic contig and genomic position                  #
# editted by Lotte van Boheemen                                         #
# last edit 4-12-14                                                     #
#########################################################################

#!/usr/bin/perl

use diagnostics;
use warnings;
use strict;

my $in = $ARGV[0];
my $pop = $ARGV[2]; #tab sep pop \t ind
my $out = $ARGV[1];
my %pop;
my %snp_to_dig;
$snp_to_dig{"A"}='1	1'; #so this is calling the keys and values in hash %snp_to_dig. Key A has value '1	1' etc.
$snp_to_dig{"T"}='2	2';
$snp_to_dig{"C"}='3	3';
$snp_to_dig{"G"}='4	4';
$snp_to_dig{"N"}='-9	-9';
$snp_to_dig{"K"}='2	4';
$snp_to_dig{"R"}='1	4';
$snp_to_dig{"W"}='1	2';
$snp_to_dig{"M"}='1	3';
$snp_to_dig{"S"}='3	4';
$snp_to_dig{"Y"}='2	3';
$snp_to_dig{"-"}='-9	-9';

my $c=1;
my %h;
my %samples;
my @samples;
my @loc;
my %label;
my $cn=1;
open POP, $pop; #this will open the population file (ARGV[2])
while (<POP>){
chomp;
my @a = split;	#this will split each line of the pop file and turns it into an array @a, $a[0] is first column in pop file, $a[1] is 2nd
$pop{$cn}=$a[1]; #calls key "1" in %pop the same as 2nd column of pop file, in this case individual
$label{$cn}=$a[0]; #calls key "1" in %label the same as 1st column of pop file
$cn++; #adds 1 to cn, so calls key 1,2,3 etc
}

close POP;

foreach my $i (keys %pop){
print "$i $pop{$i} \n";
}

open IN, $in;

#foreach (keys %snp_to_dig){
#	print "keys $_ $snp_to_dig{$_}\n";
#} # this would print off all the keys with their values in %snp_to_dig

while (<IN>){
chomp;
my @a = split; #@a is the infile, with each item seperated by the command 'split'
if(/^#/){ #skip any header, in this case CHROM POS <sampleIDs>
next;
}
my @tmp=();
my $loc=shift (@a); #$loc is CHROM
my $loc2=shift(@a); #$loc2 is POS
push(@tmp, $loc, $loc2); #@tmp holds name with genomic contig and genomic position
my $tmp=join '__', @tmp; #$tmp now holds CHROM__POS
push(@loc, $tmp);
my $ind=1;
foreach my $i (@a){#@a on this point only holds genotype info
print "$snp_to_dig{$i}\n";
push(@{$samples{$ind}}, $snp_to_dig{$i});
#the value of each key (ind, 1 to n) in %samples is $snp_to_dig{$i}.
# $snp_to_dig{i} is each  genotype encountered in the infile.
#So looping through @a, when encountering for example  A, this part will be '1	1'.
$ind=$ind+1;
}

}
close IN;
open OUT, ">$out";

print  OUT "label\tpopulation\t";

foreach (@loc){#print all the previous called CHROM__POS to the outfile
print OUT "$_\t"; #print out each locus after label\tpopulation (seperated by tabs)
}
print OUT "\n";#and go to next line

#this assumes each individual is sorted in the same way in the input file as in the popfile.
foreach my $s (sort keys %samples){
if ($pop){
if ($pop{$s}){
print OUT "$label{$s}\t$pop{$s}\t";
foreach (@{$samples{$s}}){
print "pop $s\n";
print OUT  "$_\t";
}
print OUT "\n";
}
}else{
print OUT "$label{$s}\t";
foreach my $t (@{$samples{$s}}){
print OUT "$t\t";
}
print OUT "\n";
}
}
```

In command line:
``` perl ~/scripts/SNPtable2structure_combined_lotte.pl snptableUG.tab.192.table.contig.random rand.structure popind```

It prints out the rand.structure table with per line the following:
popID \t sampleID \t \t Chrom_pos1 \t Chrom_pos2 etc.


I am not happy with this file as somehow I lost the 2nd allele. Also, I am pretty sure that the script assumes that the individuals in the popind file are sorted in the same way as the header of snptableUG.tab.192.table.contig.random.

```
cat rand.structure | wc -l 
``` shows that it didn't procude 2 lines per individual...

Looking at Ander's R script (structure_input_file.R), I should be able to melt and cast the table into the right format. For this, I would first melt each individual to the rows (and proceed as in your example). I've never used these functions before so not quite sure how to do this for all the 384 individuals. Also, I would not know how to link the populationIDs and translate the genotypes to digital alleles.

I can't yet get it to work without providing pop information, but I reckon that IF you want to run this into STRUCTURE without any prior population knowledge, just say all pops are 1.
Or something similar. 

 
structure -m mainparams

No admixtrure model: assign individuals/units
Admixture model: assigns genotypes

if alpha = 1 = flat prior

Run replica runs for every chosen K, at least 20 per K.

Output from STRUCTURE cannot be used in R, needs to be parsed by using the python script parse_struc.py

Use long burning times (min 100,000)

Evanno etal 2005: Detecting the number of clusters in Structure, check how fast the log-likelihood changes from K to the next K.