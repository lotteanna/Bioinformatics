**STRUCTURE**
===

http://pritchardlab.stanford.edu/structure.html

**A)    Prepare a table with randomly selected unlinked SNPs for STRUCTURE**

Previously I've worked with a script that randomely selects one locus per contig. However, some contigs might be small, and we don't want those to have an equal chance of being represented in our database. It would be better to first find SNPs at random (need to specify a maximum), and then make sure they are unlinked by only selecting unique contigs. For this we will use the file with most stringent filtering (50% of more missing data per SNP is filtered out), as we are only interested in population structure. Both STRUCTURE and Bayenv aren't able to '


Below command will shuffle the table and pick the top 2000 lines from the shuffled table (excluding the header, which starts with '#'. This is then passed on to the output.

```
module load coreutils
shuf -n 1000 snptableUG_pass.tab.p240_table | grep -v '#'  > snptableUG_pass.tab.p240.random
```
STRUCTURE can't deal with linked loci, so need to select at random one SNP per contig. 

```
perl /nfs/home/hpcsci/lotteanv/scripts/random_one_per_locus_combined_lotte.pl snptableUG_pass.tab.p240.random
```

This command isn't useful when there are many contigs, as the chances are very high that previous code already selected unique contigs (just by chance). This can be checked by counting the number of lines of both files

```
awk 'NR!=1{print $1}' snptableUG_pass.tab.p240.random | sort | uniq -c | wc -l
awk 'NR!=1{print $1}' snptableUG_pass.tab.p240.random.random | sort | uniq -c | wc -l
```

First, I lost my headers in the snp_coverage file, so I grep them from another table:

```
cat filtered_passed_snpsb.vcf | grep '#' > header
cat header snptableUG_pass.tab.p240.random.random > snp240_1000.random
```

**B)    Convert the SNPtable to a STRUCTURE readable format**

 NOTE: the order of the digital numbers is mixed up (eg AC can now be read in by Structure as CA. This means linkage mapping is impossible with this data (see Structure documentation). This is caused by the coding used initially, which is by vcf2vertical_dep_GATK_UG.pl

First, make a population document (popind) for the conversion script, which is a tab separated file. The following seems to work (don't enter the headers):
pop        sampleID
010908-1	2-010908-1-17
010908-1	2-010908-1-9
010908-1	1-010908-1-27
etc

The original file I have to translate is called snptableUG.tab.192.table.contig.random. This file is structured as follows:

Columns:
1: Chrom#
2: Locus Location
3 - 386: SamplesIDs

This table doesn't neatly have the genotypes spelled out as AA, AT, GC etc, but uses the combined coding for genotypes: A, N, W, etc.

- Digitalise and separated genotypes for each allele ('A' = 1,1; 'N' = -9,-9; 'W' = 1,2 etc)
- Concatenate Chrom and Locus Location
- Add in population information for each individual (through popind file)
- Transpose columns and put them in the right order for STRUCTURE
- Optional: Write data for each individual in 2 rows (ONEROWPERIND=0 if divided in 2 rows)



Use R script 

```
### SNPTable2STR

### History
#7/Jan/2015 - fixed bug found by Lotte where if the first individual was found to be missing
## then the script would crash. The issue was that when that happened, the outputGenos object
## was not created when i == 1 in the loop, so when the i > 1 came around, the ifelse statement
## test would lead to FALSE, and it would expect the outputGenos object to exist. The ifelse
## statment has now chaned to test if the object outputGenos does not exist. If it doesn't
## then it creates it.

#### Auxiliary functions #####

#conversion table

genotypeCodes = list("A"=c(1,1),
"T"=c(2,2),
"C"=c(3,3),
"G"=c(4,4),
"Y"=c(2,3), #C or T
"R"=c(1,4),  #A or G
"S"=c(3,4), #G or C
"W"=c(1,2), #A or T
"K"=c(2,4), #G or T
"M"=c(1,3), #A or C
"N"=c(-9,-9),
"-"=c(-9,-9))

# xref inds and pops

xref_ind = function(popinfo,indlist){
popinfo = popInfo
indlist = inds
#figure out number of populations and number of samples per pop
tallyPops = table(popinfo$pop)
nPops = length(tallyPops)

#number the populations
popIndex = 1:nPops
names(popIndex)<-names(tallyPops)

#dictionary of populations
popIndDic = popinfo$pop
names(popIndDic) = popinfo$ind

return(list(pIx=popIndex,pDic=popIndDic))
}

description = function(){
cat(paste(rep("#",80),sep="",collapse=""),"\n")
cat("usage:\n")
cat("snpTable2str.R <snptable> <poptable> <outfile>\n")
cat("\n\n")
cat("This script takes three input parameters:\n")
cat("\t\tsnptable: a filename. File contains a SNP genotype matrix (locus x individuals)\n\n")
cat("\t\tpoptable: a filename. File contains a tab delimited table mapping individuals to populations (pop,ind)\n\n")
cat("\t\toutfile: a filename. File to save structure output. Two lines per individual, plus one row of loci names (identified as CHROM_POS)\n\n")
cat(paste(rep("#",80),sep="",collapse=""),"\n\n")
}

####### main script ########

args <-commandArgs(trailingOnly=T)

if(length(args) != 3){
description()
} else {
genTable = read.table(args[[1]],header=T)

popInfo = read.table(args[[2]],header=F)
names(popInfo) = c("pop","ind")

#ind index x-ref
## individual genotype list
inds=gsub(pattern = "X",replacement = "",names(genTable)[3:length(names(genTable))])
inds = gsub("\\.","-",inds)

xrefIndex = xref_ind(popIndex,inds)

outfile = args[[3]]

#create locus id line
loci = paste(genTable$CHROM,genTable$POS,sep="_",collapse="\t")
cat(loci,"\n",file=outfile)

#create a vector to catch individuals that fail to cross reference to a population
missingPop = c()

#counter for total number of individuals successfully converted and saved to file
nInd = 0

#loop through columsn of genTable, create one big string with all converted genotypes
# this will save in time by saving to file only once
for(i in 1:length(inds)){
recodedGenotypes = matrix(unlist(genotypeCodes[as.character(genTable[,i+2])]),nrow=2)
ind = inds[i]
pop = as.character(xrefIndex$pIx[xrefIndex$pDic[inds[i]]])

#check if pop cross-reference failed
if(is.na(pop)){
missingPop = c(missingPop,ind)
next
}

nInd = nInd + 1

strGeno=paste0(paste(paste(ind,pop,sep="\t"),apply(recodedGenotypes,1,paste,collapse="\t")),collapse="\n")
outputGenos=ifelse(!exists("outputGenos"),strGeno,paste(outputGenos,strGeno,sep="\n"))
}

#save genotypes to file
cat(outputGenos,"\n",file=outfile,append=T)

#finish things off
cat("Conversion finished\n")
cat("Converted ",as.character(nrow(genTable))," loci and ",as.character(nInd)," individuals\n")
if(length(missingPop)>0){
cat("A total of ",length(missingPop)," individuals had no pop assignment and were skipped:\n")
for(i in missingPop){
cat("\t\t",i,"\n")
}
cat("The list is saved in the file missingPop.txt\n")
write.table(x = missingPop,file = "missingPop.txt",row.names = F,col.names = F)
}
cat("\n")
}

```

Run this by typing (make sure to have R version 2.15 or higher loaded):

```
module load R/2.15.0
Rscript snpTable2str.R snptableUG_pass.tab.48.random.random ../popind ../test_infile.stru
```

---

**C) Setting up parameter files**

Now, we need a set of mainparams and extraparams to start some runs. These parameters depend on the input file, as well as several assumptions made by the user. These assumptions depend on the situation (i.e. admixture, no admixture) or require tweaking by the user.

Helpful information:
Evanno et al 2005: Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study.

Porras-Hurtado et al 2013: An overview of STRUCTURE: applications, parameter settings, and supporting software

http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/structure_doc.pdf

>   No admixture model: assign individuals/units
>   Admixture model: assigns genotypes
>   if alpha = 1 = flat prior
>   Use long burning times (min 100,000)
>   Run replica runs for every chosen K, at least 20 per K

---

**D)	Running STRUCTURE**
Great! Now start some runs, preferably in parallel. Ideally, we want to work our way up the ultimate number of iterations and long burnin times, but start with some subsamples and small runs to evaluate the data.

Option 1:
```
#!/bin/sh

module load structure/2.3.4

START=$(date +%s)

maxRep=$1

for i in $(eval echo "{1..$2}")
do
    task_id=$i
    echo "This is task number $task_id\n"
    K=$[1+($task_id/$maxRep)]`
    rep=$[($K*$maxRep)-$task_id]
    out=`echo output_K"$K"_r"$rep"`
    chain=`echo chain_K"$K"_r"$rep"`
    seed=`eval od -vAn -N4 -tu4 < /dev/urandom`
    structure -K $K -o $out -D $seed > $chain
    END=$(date +%s)
    DIFF=$(echo "$END - $START" | bc)
    echo "Finished job K=$K, replicate=$rep, and took $DIFF seconds. It had seed $seed" >> log
done

```

Run this on the command line using
bash par_struc.sh <max number of repetitions> <max number of tasks>

OR Option 2 (much better option)

log a job to run it parallel:

```
#!/bin/sh
#$ -S /bin/sh
#$ -l h_rt=24:00:00
#$ -cwd
#$ -M XXXX@monash.edu
#$ -m beas
#$ -cwd
task_id=$[$SGE_TASK_ID-1]

K=$[1+($task_id/$maxRep)]

rep=$[($K*$maxRep)-$task_id]

out=`echo output_K"$K"_r"$rep"`

seed=`eval od -vAn -N4 -tu4 < /dev/urandom`

echo $task_id
echo $total_tasks
echo $K
echo $rep
echo $out
echo $seed

module load structure/2.3.4

structure -K $K -o $out -D $seed
```

qsub -v maxRep=20 -V -t 1-100 struc.job

---

**E)    Make output readable for R**

Output from STRUCTURE cannot be used in R, needs to be parsed.

```
#!/bin/sh
search=$1
out=chains_$search.txt
files=`find . -name "struc_au.o*" -exec grep -l "$search" {} \;`
n=1
for f in $files
do
if [ "$n" -eq "1" ]
then
grep -m1 "Rep" $f | sed 's/#://g ; s/[\s^\S\n]/\t/g' >$out
fi
sed -n '/BURNIN/,/MCMC complete/{ s/://g; s/[\s^\S\n]/\t/g; /^[0-9]/p}' <$f >>$out
n=$(($n + 1))
done
```

Run with:

```
bash parse.sh K2_
```

Note: be careful when parsing K1, it will also parse K2, so better type K1_ and change 
name of output file

This will produce chains_K2_.txt etc

The output gets really big if the number of runs and iterations are large. I have 
implemented my 'manual thinning', so I can analyse the runs quicker. All this is doing is
 removing datapoints, so from 1M points I go to 1,000, which is still a fair amount to check data.


```
for i in chains_*.txt;
do awk '!(NR % 1000)' $i > ${i/.txt/}1000sum.txt; done
```


**F)    Run structureHarvester on structure outputs**

http://users.soe.ucsc.edu/~dearl/software/structureHarvester/


The output folder should also hold files that look similar to: ```output_K8_r6_f```. These are the files structureHarvester needs. Apparently this program is a little picky, and needs the input dir to be 'clean'. Move all the files that are direct outputs from the struc.jobs to job_files:

```
mv struc_au.* job_files
```

Make sure to also place harvesterCore.py in the same directory as structureHarvester.py
NOTE: this script has to be WITHIN the dir where the in and output dirs are located, otherwise it won't run. 
Make sure to load python version 2.6.6

```
module load python/2.6.6 
python structureHarvester.py --dir=structure_output --out=clumpp_files --evanno --clumpp
```


OR

Zip all the *_f files (produced by STRUCTURE) and upload to
http://taylor0.biology.ucla.edu/structureHarvester/


**G)    Check convergence of chains**

This output of pstruc_parse.sh and structureHarverster.py can be checked with structure_analysis.R

```{r}
#structure run analysis
library(ggplot2)

#first run the parse_struc.py script

#load chain
chain_k2 = read.table("chain_K2_sum.txt",header=T,na.strings='-')

#plot alpha chains across Reps
ggplot(chain_k2,aes(x=Step,y=Alpha,col=Rep))+geom_line()

#check alpha histograms
ggplot(chain_k2,aes(x=Alpha))+geom_histogram()+facet_grid(Rep~.)

#plot LnLike chain
ggplot(chain_k2,aes(x=Step,y=Ln_Like,col=Rep))+geom_line()

#check Ln_Like histograms
ggplot(chain_k2,aes(x=Ln_Like))+geom_histogram()+facet_grid(Rep~.)

#check correlation between F
ggplot(chain_k2,aes(x=F1,y=F2))+geom_point()+coord_fixed(0.5)+facet_grid(Rep~.)

################################################################
# plotting evanno data after applying structureHarvester.py

evanno_res = read.table("evanno.txt",header=F,comment.char='#')
names(evanno_res) = c("K","reps","mean_LnPK",	"sd_LnPK",	"Ln1K",	"Ln2K",	"Delta_K")

ggplot(evanno_res, aes(x=K, y=mean_LnPK)) + 
		geom_errorbar(aes(ymin=mean_LnPK-sd_LnPK, ymax=mean_LnPK+sd_LnPK), width=.1) + 
		geom_line() + geom_point()

ggplot(evanno_res, aes(x=K, y=Delta_K)) + geom_line() + geom_point()

```

**H) Running CLUMPP**
https://web.stanford.edu/group/rosenberglab/clumpp.html

CLUMPP expects a parameter file which is based on the number of clusters (K), the number 
of individuals or populations – depending on you're running datatype 0 or 1, 0 being for
the indfiles and 1 for the population files – the number of runs (I think determined by 
the number of runs you ran in Structure) and other options. This parameter file will call 
output files produced by structureHarvester, K1.indfile, K1.popfile etc.


In command line type
```
module load clumpp
CLUMPP param_ind
CLUMPP param_pop
```

**I)Running Distruct**

Rename the .output files from CLUMPP to .indivq and .popq for the indivual and population 
run respectively

```

```

ps K2.ps

By ordering the labels in K*.names and K*.languages, you can change the order in the plot.
By setting PRINT_INDIVS to 1, you print the q-scores per individual. Setting this to 0 will print the q-scores per population
In K*perm, it is possible to specify the colours. The first colour in the list will be 
drawn first in the bar (bottom). 
The numbers in front will specify something with the q-value. 
It is really just messing around with these numbers to get the desired order of colours








---
Parse STRUCTURE output chains in python (not working)

When using *option1* in *D*, use the python script parse_struc.py

```
#!/usr/bin/env python
'''
parse_struc.py <chain_prefix>

e.g.,

./parse_struc.py chain_K2
'''

import sys,os,re,string

prefix = sys.argv[1]

def find_files(prefix):
'''
will return a list of files in the current directory that contain the prefix
'''
return [f for f in os.listdir('.') if re.search(prefix,f)]

def parse_files(filelist,prefix):
'''
output a file prefix_sum.txt with the summary chains
'''
fo = open(prefix+'_sum.txt','w')
file_count = 0
out = []
for f in filelist:
print f
fi = open(f,'r')
rep = f.split('_')[-1]
tmp1 = []
count = 0
start = 0
mycounter = 0
for line in fi:
tmp = line.strip().lstrip()
if re.search('starting MCMC',tmp) and start == 0:
start = 1
if re.search('[0-9]{3,}:',tmp) and start == 1:
tmp1.append(re.sub(':','',tmp))
else:
if re.search('BURNIN completed',tmp):
count = 1
else:
if count==1 and re.search('Rep#:',tmp) and file_count==0:
header=tmp
header = re.sub('[ ]{2,}',';',header)
header = re.sub('[ ,]{1}','_',header)
header = re.sub(';',' ',header)
header = re.sub('Rep#:','Step',header)
header = 'Rep '+ header
count +=1
file_count+=1
else:
continue
#tmp1 = [re.sub(':','',r.strip().lstrip()) for r in fi if re.search('[0-9]{3,}:',r)]
tmp2 = [re.sub('[ ]{1,}',' ',r) for r in tmp1]
tmp3 = [re.sub('--','- -',r) for r in tmp2]
tmp4 = [rep+' '+r+'\n' for r in tmp3]
out.extend(tmp4)

fo.write(header+'\n')
for r in out:
fo.write(r)
fo.close()
return out

if __name__=="__main__":
flist = find_files(prefix)
parse_files(flist,prefix)
```

Now in the directory where all the outputs are:

```
python ~/scripts/parse_struc.py chain_K1_ #(etc, for every K)
```

This will create a *sum.txt file for every K







-------
Run the conversion script in Perl (not working)
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
```
perl ~/scripts/SNPtable2structure_combined_lotte.pl snptableUG.tab.192.table.contig.random rand.structure popind
```

It prints out the rand.structure table with per line the following:
popID \t sampleID \t \t Chrom_pos1 \t Chrom_pos2 etc.

I am pretty sure that the script assumes that the individuals in the popind file are sorted in the same way as the header of snptableUG.tab.192.table.contig.random.
Also, this script prints out


```
cat rand.structure | wc -l 
```

shows that it didn't produce 2 lines per individual...


-------------------



