Filtering the vcf files
===

Useful links

```
http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
https://www.broadinstitute.org/gatk/guide/article?id=3225
```

The vcf file is beautiful, but contains some info we're not that interested in, e.g. bad mapping qualities, loci with very low or high coverage, extreme heterozygosity, lot's of missing data, alleles with very low frequency etc. So we're going to filter that all out.


**A)Retrieve only raw SNPs from vcf file**

I don’t spend too much time on indels, as I am more interested in SNPs and it is uncertain in this case how reliable indels are. Below steps however, could also be followed for retrieving and filtering indels. This requires different treshholds, see Broad Institute website.

```
#gatk_selectSNP.sh
module load java/jdk1.7.0_21
module load gatk/3.1.1

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ~/ragweed/GBS/raw_common/ref/ragweed_10Feb2016_2ABsE_uppercase70.fasta.fa \
-V snps.raw.vcf \
-selectType SNP \
-o raw_snps.vcf
```

NOTE: I will still use grep -v 'INDELS', which will get everything BUT the indels, to make absolutely sure I don't have any indels


**B)Filter vcf file**

```
# vcf_gatk_hardfilter.sh
module load java/jdk1.7.0_21
module load gatk/3.1.1

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/ref/ragweed_10Feb2016_2ABsE_uppercase70.fasta.fa \
-V snps.raw.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || MQ < 40.0 || ReadPosRankSum < -8.0" \
--filterName "my_snp_filter"  -o filtered_snps.vcf
```

**C)Futher filterin vcf file. Custom edition**

Now we will locate called SNPs and genotypes and filters on set quality and depth it will produce an 'N' for every position (SNP/individual) for which no genotype was
called by UG or failing the filter we set.

NOTE: you are now entering in the domain of custom scripts. Although you will be working with your data more hands-on, you will also produce tables that are custom-made. This means that there are not necessarily easy format scripts to manipulate the data (as you need to do this for every-single-program)

The script used is from the collection of Kay Hodgins, and the description of set thresholds originates from her blog at http://www.zoology.ubc.ca/~rieseberg/RiesebergResources/?p=25021

<cite>
qual: This is the phred-scaled probability of a SNP occurring at this site. A score of 20 means that there is a 1 in 100 chance that the SNP is a false positive. A score of 30 means that there is a 1 in 1000 chance of a false positive. There is a Qual score for each variant site.

MQ: This is the phred-scaled probability that the read is mapped to the correct location (a low map score will occur if reads are not mapped uniquely at that site– i.e. they come from a region that is repeated in the genome). There is a MQ score for each variant site.

GQ: This is the phred-scaled probability that the genotype being called is correct, given that there is a SNP at that site. There is a GQ for each individual.

Minimum and maximum individual read depth: Sometimes I have found genotypes being called based on a small number of reads (e.g. 5) although the GQ is relatively high (>20). Therefore, I will likely increase the GQ threshold or also have a minimum depth requirement. Also, high depth indicates that there are repetitive regions aligning to that site so the SNP may not be real.
</cite>


```
cat filtered_snps.vcf | grep -v 'INDEL' | /nfs/home/hpcsci/lotteanv/scripts/vcf2vertical_dep_GATK-UG.pl > snptableUG.tab
```

Get the number of filtered SNPs:
```
cat snptableUG.tab | grep Scaf* | wc -l
```

vcfdepth_lotte.pl is the exact same as vcf2vertical_dep_GATK-UG.pl but will produce
'depth.txt', which gives the depth of the # of reads for each called genotype

Explore the depth distribution with R
```
R
d<-read.table('depth.txt')
head(d) #R called my column V1
hist(d$V1)
```

**D)Filter SNPtable on population genetic parameters**

Get summary statistics and filter SNP table based on minor allele frequency, heterozygosity and missing data using snp_coverage.pl

```
perl /nfs/home/hpcsci/lotteanv/scripts/snp_coverage.pl snptableUG.tab 
```

Minor allele frequency: Low frequency SNPs could be due to errors and are not useful for
outlier tests and several other tests of selection (although they are for site frequency
spectrum tests)
Heterozygosity: High or fixed heterozygosity could indicate parology.
Missing data: SNPs that have low coverage are removed
explore data in R using the summary output:

```
d<-read.table('snptableUG.tab.summary',col.name=c("scaf_no","position","allele_tot","genotype_tot","top_o","sorted","hash_0","mj","sort_1","hash_1","mn","het"))
hist(d$genotype_tot,xlim=c(0,50),breaks=1000)
hist(d$het)
hist(d$mn) 
```

Once decided on cut-off values (I used mn_cut = 0.05, het < 0.7 and geno_cut =
<half of total number of individuals>)
the snptableUG.tab.table is the new filtered output file

Check quality per individual
add header names to snptableUG.r.tab.table by deleting # from header in snptableUG.tab to
new file

```
R
d<-read.table("snptable.tab.r192.table", strip.white=T,header=T)
head(d)
lst<-colnames(d[,-(1:2)]) # discards the first 2 columns in the lst headers
```

Get rid of all the X in the headers and replace the . in the sampleIDs with dashes

```
lst<-gsub(pattern = "X",replacement = "",lst) 
lst<-gsub("\\\\.","-",lst) #doesn't work
```

Count the total number of N's for each individual

```
countN<-sapply(lst,FUN=function(x,df){sum(df[,x]=="N",na.rm=T)},df)
```

Calculate the percentage of N per SampleID divided by total number of SNPs

```
percN<-sapply(lst,FUN=function(x,d){sum(d[,x]=="N",na.rm=T)}/nrow(d),d)
```

And make a histogram

```
hist(percN)
```
write.table(percN,"percN_filt192i.txt",quote=F,sep="\t")
#write.table(countN,"countN_filt192i.txt",quote=F,sep="\t")
```
lst<-gsub(pattern = "X",replacement = "",colnames(d[,-(1:2)]))

Copy these files to local computer --> open new terminal but DO NOT ssh into MCC


Check the number of SNPs gained through each method

```
cat snptableUG.tab.table | wc -l
```
---

**OPTIONAL** only necessary when scaffolds were seperated by AAAA to reduce genome size

Put scaffolds back into original contig positions
```
perl /nfs/home/hpcsci/lotteanv/scripts/scaffold2contig.v3.pl snptableUG.tab.table ~/ragweed/WGS/soap_assembly/scaffold_order.soaprunk61.contig\
```

Check how this all went
This will get every SNP for each unique Scaffold
```
-cat snptableUG.tab.table | grep Scaffold* | cut -c 1-13 | sort | uniq -c | wc -l
```

353; this means that 35% of the arbitrary scaffolds made contain SNPs

After putting back to original contigs

```
awk 'NR!=1\{print $1\}' snptableUG.tab.table.contig | sort | uniq -c | wc -l
```
NOTE: this is wrong coding! Printing this will give the unique number of Chrom, whereas I need a combination of Chrom and Pos

Scaffolds were arbitrary, so it makes sense that the SNPs are on more contigs than scaffolds.