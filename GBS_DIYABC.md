DIYABC
===

> Code by Lotte van Boheemen

Convert vcf file to DIYABC input file via python script (note that custom filtering is skipped this way)

https://github.com/loire/vcf2DIYABC.py

Download DIYABC + manual (registration required for latter)

http://www1.montpellier.inra.fr/CBGP/diyabc/

**A)	Create SNP datafile**

DIYABC requires a datafile specific for the program. There is a python script which 
re-formats vcf files to the required SNP file. The only problem with this option is 
that I am not able to do my final step of filtering, which is taking out all the SNPs
which have a SNP call rate of ...% or less. Maybe it is possible to do this step afterwards

User needs to provide a popfile.txt which specify individuals sex and population of origin
popfile must be formated as follows:
indiv sex pop

Sex is either M or F, or 9 if undefined
```
module load python
python vcf2diyabc.py infile.vcf popfile.txt sexratio
```

So I either have to subset my vcf file, or I have to subset and rename pops in the abc 
file. Possibly the first option is better as the header might be written based on
information found in the vcf file. Then I need to check if the vcf is holding any information I don't want


**B)	Filter SNP datafile**

As this DIYABC SNP table is derived directly from a (filtered) vcf-file, custom filtering
based on SNP call rate (the percentage of individuals with called genotypes at a given
SNP) was skipped. Still, this is an important step, so I made a short R script to filter
the DIYABC file. For explanations on filtering, see GBS5_vcffilter.md

In the SNP datafile, 0 = homozygous genotype for the non reference allele, 1 = heterozygous
genotype for the reference allele, 2 = homozygous genotype for the reference allele. We
thus need to count the number of '1' and delete all columns with a frequency of '1' > 0.7

A SNP is monomorphic when there is only 1 allelic state. Although this should not influence
any analysis, monomorphic SNPs increase computing time. I will filter out all SNPs which
are completely homozygous for one of the alleles.

This conversion script requires as input a diyabc snp file. The goal of the current script is to 
- filter out SNPs with a heterozygosity higher than a set threshold
- filter out SNPs with a set amount of missing data
- filter out monomorphic sites

---

Reformat DIYABC input file
===

> Code by Lotte van Boheemen


This conversion script requires as input a diyabc snp file. The goal of the current script is to 
- filter out loci with a set amount of missing data

Read data

```{r}
data<-read.table("filt_completeMAF.05.DIYABC.snp", header=T)
```

```{r}
geno_cut <- 0.5 # SNP call rate
het_cut <- 0.7 # maximum heterozygosity

# create boolean values (TRUE/FALSE) for each column indicating if they passed the filter
het_pass <- apply(data, 2, function(x) length(which(x == "1"))/length(x)) <= het_cut
null_pass <- apply(data, 2, function(x) length(which(x == "9"))/length(x)) <= geno_cut
homref_pass <- apply(data,2, function(x) length(which(x == "0"))/length(x)) != 1
homalt_pass <- apply(data,2, function(x) length(which(x == "2"))/length(x)) != 1

#check if this worked by counting the number of loci which passed each filter
length(which(het_pass== "FALSE"))
length(which(het_pass== "TRUE"))
length(which(null_pass== "FALSE"))
length(which(null_pass== "TRUE"))
length(which(homref_pass== "FALSE"))
length(which(homalt_pass== "FALSE"))


# write columns passing the filter to new data frame
#trydata <- data[het_pass]
trydata2 <- data[null_pass] #this works
#sdata <- data[null_pass | het_pass] #but this doesn't

# Changing the order around of filtering
het_pass <- apply(trydata2, 2, function(x) length(which(x == "1"))/length(x)) <= het_cut
#null_pass <- apply(trydata, 2, function(x) length(which(x == "9"))/length(x)) < geno_cut #for some reason I'm getting those with too much missing data
#sdata <- data[null_pass]
sdata <- trydata2[het_pass=="TRUE"]
ncol(trydata2)
ncol(sdata)

# Maybe it's a good idea to sort on MAF before putting through diyABC. If I am going to subset before, MAF will be different. 

# write data drame to file
write.table(sdata,"diyabc_snp50MAF.05_all")
```


Subsetting a random number of columns as follows:
```{r}
subset<-cbind(sdata[,1:3],sdata[,sample(ncol(sdata), 1000)])
write.table(subset,"diyabc1000_snp50MAF.05_all")
```

This will give every column a number and we don't need this. So find all characters 
starting with "A." and replace with A:

```
sed -E ’s/A.[^ ]*/A/g’ diyabc_snp50MAF.05_all > test
```

Remove all quotes in the file

```
sed 's/\"//g' test > diyabc_snp50MAF.05_all
```

*C)	Recode populations**

To run ABC, we need to look at very simplistic scenarios. DIYABC is not happy about being
given input which it won't use, so we need to make all the subsets based on recoded
population identifiers before loading data into the program. This way, we can also remove
any populations which we are not going to analyse (e.g. outliers).


At this moment, I edit the populations in excel to add/remove populations or assign 
clusters according to STRUCTURE. This can be done in R as well, but I haven't gotten
around to scripting this yet.

To get the best format, copy out of excel into TextWrangler and save as .snp

