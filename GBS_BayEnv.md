BAYENV ===


**This script requires climate data. Produce this data with
QGIS_coordinate_climatevalue_Walkthrough.md**

BAYENV:
1) Variance covariance matrix with subset, multiple runs  
2) Check convergence MCMC runs/correlation matrices (cov2cor in R)  
3) XTX on subset —> identify outliers  
4) remove outliers original subset and repeat 1&2  
5) Env/xtx  on complete dataset (using corrected variance/covariance matrix)  


A) **Transform SNPtable to Bayenv format**

To transform the custom-filtered SNPtable to Bayenv(2) format, we need a file showing the
individual, followed by population. Preferably in the same order as the SNPtable. It can look
something like:

2-010908-1-17   010908-1 
2-010908-1-9    010908-1 
2-160808-1-11   160808-1 
2-160808-3-11   160808-3

Then convert with the perl script (authorship unknown, from shared Hodgins drive

```
#SNPtable2bayenv_memory_lotte.pl
#!/usr/bin/perl

use warnings;
use strict;

unless (@ARGV == 3) {die;}

my $in = $ARGV[0]; #merged bsnp table loci (rows) pop (columns)
my $pop = $ARGV[1];# list of population with ind \t pop
my $out = $ARGV[2]; #name of out file
my @pop=();
my $loc;
my %pop;
my $c=1;
my %h;
my %counter;
my %samples;
my @samples;
my @loc;
my $locNum;
my @a;
my $pop_count;
my %pop_no;
my @ind;
my $cn=1;
my $num=1;
open POP, $pop;
while (<POP>){
        chomp;
        unless (/^\s*$/ | /^ID/){
                @a = split;
                $pop{$cn}=$a[1];
                $pop_no{$a[1]}=1;
                $cn++;
}
}
close POP;
open OUT3, ">$out.pop_order";
my $no=1;
foreach my $p (sort keys %pop_no){
                print OUT3 "$no $p\n";
                push(@pop, $p);
                $no++;
}

foreach my $p (sort keys %pop){
                #print "pop $p $pop{$p}\n";
}

my $pop_no=scalar keys %pop_no;
print "there are $pop_no populations\n";


open IN, $in;
open OUT, ">$out";
open OUT2, ">$out.loci";

#read in snp table
while (<IN>){
        chomp; #delete hard returns
        if(/^#/){ #skip any header
                next;
        }
        @a = split;
        $locNum++;
        my @l=();
        my $l=shift (@a);
#       shift(@a);
#       shift(@a);
        my $l2=shift(@a);
#       my $l3=shift(@a);
        push(@l, $l, $l2);
    $loc = join'__', @l; #name with have transcript name, genomic contig and genomic position
        $pop_count=1;
        foreach my $i (@a){
                if($i =~ /A|T|G|C/ ){
#                       print "population $pop{$pop_count}\n";

                        $counter{$loc}{$pop{$pop_count}}{$i}++; #hash of hash of a hash locus -> pop -> allele count
                        $counter{$loc}{$pop{$pop_count}}{$i}++;
                }if($i =~ /K|W|Y/){
                        $counter{$loc}{$pop{$pop_count}}{'T'}++;
                }if($i =~ /K|R|S/){
                        $counter{$loc}{$pop{$pop_count}}{'G'}++;
                }if($i =~ /R|W|M/){
                        $counter{$loc}{$pop{$pop_count}}{'A'}++;
                }if ($i =~ /M|S|Y/){
                        $counter{$loc}{$pop{$pop_count}}{'C'}++;
                }
                $pop_count=$pop_count+1;
        }
                foreach my $k (keys %counter){
                                my %alleles=();
                                for (my $i=0; $i < @pop  ; $i++){
                                                foreach my $k2 (keys %{$counter{$k}{$pop[$i]}}){
                                                                #print "$k2\n";
                                                                $alleles{$k2}=1;#get the alleles for each loci
                                        }
                                }
                                foreach my $key (keys %alleles){
                                        #       print "alleles $key\n";
                                }


                for (my $j=0; $j < @pop; $j++){
                        foreach my $key (keys %alleles){
                                unless(exists $counter{$k}{$pop[$j]}{$key}){
                                        $counter{$k}{$pop[$j]}{$key}=0;#insert 0 for each population where that allele does not exist
                                        #print "MISS $k $key $counter{$k}{$pop[$j]}{$key}\n";
                                }
                        }
                }

        }

                my @temp1=();
                my $allele_count=();
                my @temp=();
                #print out allele count and number in each population
                foreach my $key (sort keys %counter){
#               print OUT "\n";
                                $num++;
                                @temp1=();
                                $allele_count=();
                                @temp=();
                                foreach my $key2 (sort keys %{$counter{$key}}){
                #                               print "pop counter $key2 \t";
                                                $c=1;
                                                foreach my $key3 (sort keys %{$counter{$key}{$key2}}){
                                                #               print "counter $counter{$key}{$key2}{$key3}\t allele key $key3 pop $key2\n";    
                                                                $allele_count=scalar keys %{$counter{$key}{$key2}};
                                                                if($c==1){
                                                                                push(@temp, $counter{$key}{$key2}{$key3});
                                                                }if($c==2){
                                                                                push(@temp1, $counter{$key}{$key2}{$key3});
                                                                }
                                                                        $c++;
                                                }
                                }
  }
                                }
                                if($allele_count == 2){#do not print if there are more than 2 alleles
                                                print OUT2 "$num\t$key\n";#print out the loci order for those that have 2 alleles
                                                foreach (@temp){
                                                print OUT "$_\t";
                                                }
                                                print OUT "\n";
                                                foreach (@temp1){
                                                                print OUT "$_\t";
                                                }
                                                print OUT "\n";
                                }


                }
%counter=();
}


```

Call script SNPtable2bayenv_memory_lotte.pl

``` perl ~/scripts/SNPtable2bayenv_memory_lotte.pl snp240_1000.random indpop snp240_bayenv ```


This might not work, possibly in the case that the snptable has more individuals then the indpop
file. In the SNPtable2structure reformatting script this was not an issue, as the script also output
a list of "missingPops.txt". There is of course the option to modify the SNPtable2bayenv script to
allow for this discrepancy, but my choice was to use PGDspider to change the input file for
Structure. http://www.cmpg.unibe.ch/software/PGDSpider/. One thing to note though is that you might
want different input files for Structure and Bayenv. For Structure, I applied very stringent
filtering, with <50% missing data per locus. For Bayenv though I allowed for missing data up to 90%
per locus, as I need a dataset more representative of the entire dataset (cite, Coop?). Make sure to
put in the correct format for Structure in PGDspider. The script in GBS6_STRUCTURE will output a
file with diploid data on 2 consecutive rows, no phase information, missing value code of -9, SNP
format, inclusion of locus, individual names and population identifier. There is no additional
information in the file.


**B)    Estimate covariance matrix**

First, we need to make the covariance matrix. This will take a long time.

In below command (for BAYENV), the first '0' indicates we want to estimate the covariance matrix,
'89' is the number of populations, '-23457' is the negative seed, '10000' is the number of
iterations (100,000 is standard for covariance matrix), 'snp48_bayenv' is the snpstable we want to
estimate the covariance matrix from. Note this is a random subset of non-linked SNPs, as described
in the structure walkthrough. According to the Bayenv2 manual, increasing the number of SNPs (e.g.
thousands) doesn't result in large differences in estimation of the covariance matrix. This SNPtable
is different from the SNP table we will use to do further calculations in bayenv, which will be more
inclusive.

``` bayenv 0 89 -23457 100000 snp48_bayenv > matrixall.out ```

For BAYENV2, the command is a little different (-i is input file, -p is number of populations, -k is
number of iterations, -r is the root, this can be any number). This will print to "matrices" the
updated covariance matrix for every 5000 iterations:

``` bayenv2 -i SNPSFILE -p 89 -k 100000 -r 63479 > matrices ```

Either can be submitted as an array job. Note that the computing time specified in the job is valid
for every array job, not for all array jobs together. If you want to submit 100 array jobs, the
command is

``` qsub -t 1-100 bayenv2.job ```

**C)    Run Bayenv with env file**

This step will require 3 input files, ENV, SNPFILE and MATRIX_FILE. The latter 2 look the same as
files produced/used above, but they are not! The code to run this is as follows (don't run yet, I am
going to explain every step):

``` bayenv2 -i SNPFILE -m matrix_file -e env -p 4 -k 1000 -n 8 -t -r 429 ```

-p is the number of populations, -k is the number if iterations (needs to be larger then 1000 or it
will divide by 0 in the bayenv source code), -n is the number of environmental variables, -t
indicates 'rest' mode, and will make sure the Bayes Factor is calculated for every SNP, -r is the
random root.

1]   Environmental variables are standardised by substracting the mean and dividing through the
standard deviation. Even though I will be analysing subsets for each of the continents and all of
the samples, I won't subset the environmental variables and standardize using data within each
range. If my goal is to look at parellelism between introduced ranges, I need to hold the gradients
similar. I cannot compare if standardisation is different.

The format of the environmental data for 4 populations and 8 environmental factors is shown below.
Each column is a population, in the same order as the SNPSFILE and input file.

-2.053517756	-1.905325197	-1.828770599	-1.836090653
1.563096989	1.594375406	1.600992044	1.574604643
-1.181446926	-1.171381351	-0.980135416	0.378717279
-1.465018067	-1.47330868	-1.550562123	-0.793855224
0.614591847	0.56667764	1.892304032	-0.455492108
-0.572885276	0.313220562	0.060047465	0.249927288
0.482329177	1.521168194	2.759319348	-0.571609245
0.872405436	1.162237909	2.359372036	-0.627162575


2]  The matrix_file is a summary or excerpt from above "matrices". Either average among all computed
covariance matrices (not done here, needs evaluation of convergence), or use the last one. 
Note that bayenv is VERY sensitive in the use of whitespaces vs tabs, and a simple copy 
paste might change tabs into whitespaces. To replace
multiple spaces with a single tab, use the following:

``` cat spaced-file | sed 's/ \+/\t/g' > tabbed-file ```

Another option is to yank the lines out if the matrix file. Go to the first line of the last matrix
<linenumber>gg and type the number of lines you want to copy (should be equal to number of populations in your file
<numberoflines>yy. Go to a new file and type p.    
If there are more then 50 lines, we need to adjust vim to store more in the registers. Type
```:set viminfo='20,<1000,s1000```
Followed by the yanking.

Open the new file and simply type <p>

Now, you want to average over multiple covariance matrices from diffferent runs. Paste all the
funal matrices of each run into a new file and load R to get a mean matrix:

```
module load R
A<-read.table("matrix_sel2_allr")
B<-read.table("matrix_sel_allr")
C<-read.table("matrix_sel3_allr")
D<-read.table("matrix_sel4_allr")
E<-read.table("matrix_sel5_allr")

my.list<-list(A,B,C,D,E)
meanm<-Reduce("+", my.list) / length(my.list)
write.table(meanm,"mean_matrixallr",sep="\t",col.names=FALSE,row.names=FALSE)
```



3] SNPfiles contain the allele counts (on 2 lines for diploids) for ONE SNP ONLY. This means that
SNPSFILE needs to be cut into little pieces. Below script will do this. Sidenote on this script:
authorship is unknown, copied from shared folder, edits by Simon Michnowicz

```	
#!/usr/bin/perl
my $in = $ARGV[0]; #merged bsnp table loci (rows) pop (columns)
my @array=();
print "Input file is $in\n";
open IN, $in;
while (<IN>){
        chomp;
        push (@array, $_);
#       print "$_\n";           
}
close IN;
my $len=$#array+1;
my $counter=1;
print "Size of input file is $len lines\n";

my $fileIndex;
for(my $f=0; $f < $len; $f++){
        print "Index is $f\n";
        if ($counter==1){
                $counter=2;
                $fileIndex=$f+1;
                open OUT, ">loci$fileIndex";
                my @temp=split(" ",$array[$f]);
                foreach (@temp){
                        print OUT "$_\t";
                }
                print OUT "\n";
        }elsif ($counter == 2){
                $counter=1;

                my @temp=split(" ", $array[$f]);
                foreach (@temp){
                        print OUT "$_\t";

                }
                print OUT "\n";
                close OUT;
                print "Fileanme=loci$fileIndex\n";
                my $command="bayenv2 -i loci$fileIndex -m matrix -e env -p 7 -k 10000 -n 23 -t -r 429";
                print   "$command\n";
                                system ("$command");
                }
}
#system "rm loci*";

```

Run this script as below. au42_48.bay.txt is the SNPSFILE

``` perl bayenv_split.pl au42_48.bay.txt ```

This will produce files containing information for every locus, the SNPfiles. Note that the above
script is producing the right number of files in the right order, but the naming is off, as it skips
every even number.

**D)	Running XtX matrix**

We can use the same script as above to run the XtX matrices. This also requires an environment input
file, but won't use it. The -t flag is left out, the -X flag is added. For this, we can also create a dummy environment file, filled with 0 and the amount of columns as populations

bayenv2 -i loci1a -m matrix_sel_allcra -e env -p 86 -k 100000 -n 1 -t -X -r 429

Or change in bayenv_split:
```
my $command="bayenv2 -i loci$fileIndex -m matrix -e env -p 7 -k 100000 -n 23 -t -X -r 429";
```