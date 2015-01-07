{\rtf1\ansi\ansicpg1252\cocoartf1265\cocoasubrtf210
{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural

\f0\fs24 \cf0 \
E) filter vcf file\
this script will locate called SNPs and genotypes and filters on set quality and depth\
it will produce an 'N' for every position (SNP/individual) for which no genotype was\
called by UG\
\
call vcf file and get all the SNPs (grep -v will get everything BUT the command)\
```\
cat snps.raw.vcf | grep -v 'INDEL' | /nfs/home/hpcsci/lotteanv/scripts/vcf2vertical_dep_GATK-UG.pl > snptableUG.tab\
```\
\
Get the number of filtered SNPs:\
```\
cat snptableUG.tab | grep Scaf* | wc -l\
```\
\
vcfdepth_lotte.pl is the exact same as vcf2vertical_dep_GATK-UG.pl but will produce\
'depth.txt', which gives the depth of the # of reads for each called genotype\
\
Explore the depth distribution with R\
```\
R\
d<-read.table('depth.txt')\
head(d) #R called my column V1\
hist(d$V1)\
```\
\
\
F) get summary statistics and filter SNP table based on minor allele frequency,\
heterozygosity and missing data using snp_coverage.pl\
```\
perl /nfs/home/hpcsci/lotteanv/scripts/snp_coverage.pl snptableUG.tab \
```\
\
Minor allele frequency: Low frequency SNPs could be due to errors and are not useful for\
outlier tests and several other tests of selection (although they are for site frequency\
spectrum tests)\
Heterozygosity: High or fixed heterozygosity could indicate parology.\
Missing data: SNPs that have low coverage are removed\
explore data in R using the summary output:\
\
```\
d<-read.table('snptableUG.tab.summary',col.name=c("scaf_no","position","allele_tot","genotype_tot","top_o","sorted","hash_0","mj","sort_1","hash_1","mn","het"))\
hist(d$genotype_tot,xlim=c(0,50),breaks=1000)\
hist(d$het)\
hist(d$mn) \
```\
\
Once decided on cut-off values (I used mn_cut = 0.05, het < 0.7 and geno_cut =\
<half of total number of individuals>)\
the snptableUG.tab.table is the new filtered output file\
\
Check quality per individual\
add header names to snptableUG.r.tab.table by deleting # from header in snptableUG.tab to\
new file\
\
```\
R\
d<-read.table("snptable.tab.r192.table", strip.white=T,header=T)\
head(d)\
lst<-colnames(d[,-(1:2)]) # discards the first 2 columns in the lst headers\
#countN<-sapply(lst,FUN=function(x,df)\{sum(df[,x]=="N",na.rm=T)\},df)\
#counts the total number of N\'92s for each individual\
percN<-sapply(lst,FUN=function(x,d)\{sum(d[,x]=="N",na.rm=T)\}/nrow(d),d)\
#lst<-gsub(pattern = "X",replacement = "",lst) #gets rid of all the X in the headers, don't get this to work\
#lst<-gsub("\\\\.","-",lst) #replaces the . in the sampleIDs with dashes, don't get this to work\
#calculates the percentage of N\'92s per ind divided by total number of SNPs\
hist(percN)\
write.table(percN,"percN_filt192i.txt",quote=F,sep="\\t")\
#write.table(countN,"countN_filt192i.txt",quote=F,sep="\\t")\
```\
lst<-gsub(pattern = "X",replacement = "",colnames(d[,-(1:2)]))\
\
Copy these files to local computer --> open new terminal but DO NOT ssh into MCC\
\
\
G) Put scaffolds back into original contig positions\
```\
perl /nfs/home/hpcsci/lotteanv/scripts/scaffold2contig.v3.pl snptableUG.tab.table ~/ragweed/WGS/soap_assembly/scaffold_order.soaprunk61.contig\
```\
\
Check how this all went\
This will get every SNP for each unique Scaffold\
```\
-cat snptableUG.tab.table | grep Scaffold* | cut -c 1-13 | sort | uniq -c | wc -l\
```\
\
353; this means that 35% of the arbitrary scaffolds made contain SNPs\
\
After putting back to original contigs\
\
```\
awk 'NR!=1\{print $1\}' snptableUG.tab.table.contig | sort | uniq -c | wc -l\
```\
NOTE: this is wrong coding! Printing this 
\f1\fs22 \CocoaLigature0 will give the unique number of Chrom, whereas I need a combination of Chrom and Pos
\f0\fs24 \CocoaLigature1 \
\
Scaffolds were arbitrary, so it makes sense that the SNPs are on more contigs than scaffolds.\
\
}