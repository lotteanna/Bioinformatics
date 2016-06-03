Filtering the vcf files
===

Useful links

```
http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
https://www.broadinstitute.org/gatk/guide/article?id=3225
```

The vcf file is beautiful, but contains some info we're not that interested in, e.g. bad mapping qualities, loci with very low or high coverage, extreme heterozygosity, lot's of missing data, alleles with very low frequency etc. So we're going to filter that all out.

The output of this process is 2 files filtered as:

- gatk_selectSNP.sh (remove indels) —> input is snps.raw.vcf, Output is raw_snps.vcf

- vfc_gatk_hardfilter.sh (QD < 2, MQ < 40, MQRandSum < -12.5, ReadPosRankSum < -8) —> Input is raw_snps.vcf, Output is filtered_snps.vcf

- grepping the passed SNPs —> Input is filtered_snps.vcf, Output is filtered_passed_snps.vcf

- vcf2vertical_dep_GATK-UG.pl (GTqual >=20, MQ>=20, qual >= 20, dp = 5 - 240), note that mapping quality shouldn’t be necessary here as that was already done in the previous step —> input is filtered_passed_snps.vcf, Output is snptableUG_pass.tab

- snp_coverage_p240.pl (heterozygosity =<70, MAF >= 0.05, SNPcallrate >= 50% of samples) —> Input is snptableUG_pass.tab, output is snptableUG.tab_p240.table
OR snp_coverage_p48 (heterozygosity =<70, MAF >= 0.05, SNPcallrate >= 10% of samples) —> Input is snptableUG_pass.tab, output is snptableUG.tab_p48.table


**A)    Retrieve only raw SNPs from vcf file**

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


**B)    Filter vcf file**

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

The file size should NOT have been reduced, as the only difference should be that the snps
now contain a stamp ‘Pass’ or <Filtername>

So I just got everything that has 'PASS', and make sure I get the header too. 

An excellent way is to filter the vcf file with vcftools, which has plenty of options
<link>https://vcftools.github.io/man_latest.html</link>

Need to filter vcf file with vcftools, as I cannot use the vcf2vertical script. 
In this script used the following thresholds

The description of thresholds used in this walkthrough comes from Kay Hodgins' blog at
<link>http://www.zoology.ubc.ca/~rieseberg/RiesebergResources/?p=25021</link>

<cite> qual: This is the phred-scaled probability of a SNP occurring at this site. A score
of 20 means that there is a 1 in 100 chance that the SNP is a false positive. A score of
30 means that there is a 1 in 1000 chance of a false positive. There is a Qual score for
each variant site.

MQ: This is the phred-scaled probability that the read is mapped to the correct location
(a low map score will occur if reads are not mapped uniquely at that site– i.e. they come
from a region that is repeated in the genome). There is a MQ score for each variant site.

GQ: This is the phred-scaled probability that the genotype being called is correct, given
that there is a SNP at that site. There is a GQ for each individual.

Minimum and maximum individual read depth: Sometimes I have found genotypes being called
based on a small number of reads (e.g. 5) although the GQ is relatively high (>20).
Therefore, I will likely increase the GQ threshold or also have a minimum depth
requirement. Also, high depth indicates that there are repetitive regions aligning to that
site so the SNP may not be real. </cite>


```
module purge
module load vcftools
vcftools --vcf filtered_snps.vcf --remove-indels --remove-filtered-all --minGQ 20.0 --recode --recode-INFO-all --out filtered_pass.vcf
```

--remove-indels will remove everything marked as an indel
--remove-filtered-all will remove everything not marked as PASS
--maf <float> filter for minor allele freq
--minGQ <float> filter for minimum genotype qual
--recode --recode-INFO-all --out <filename>.vcf allow for saving the output file as vcf
	while retaining all information

The info per SNP should now like like the following:

```
SC2ABSE_765.2   239455  .       G       A       811.13  PASS    AC=2;AF=0.143;AN=14;BaseQRankSum=-4.768;DP=729;Dels=0.00;FS=0.000;HaplotypeScore=1.2942;MLEAC=2;MLEAF=0.143;MQ=59.84;MQ0=0;MQRankSum=0.482;QD=21.35;ReadPosRankSum=-0.636       GT:AD:DP:GQ:PL 
```

Generate a file reporting the missingness on a per-individual basis. 
The file has the suffix ".imiss"

```
vcftools --vcf filtered_pass.vcf.recode.vcf --missing-indv --out vcftoolsPASSED
```

Generates a file reporting the missingness on a per-site basis. 
The file has the suffix ".lmiss".

```
vcftools --vcf filtered_pass.vcf.recode.vcf  --missing-site --out filtered_pass.vcf.recode.vcfPASSED
```

Additional filtering (SNP quality, min and max depth, min mapping quality)

```
cat filtered_pass.vcf.recode.vcf | vcf-annotate --filter Qual=20/MinDP=5/MaxDP=240/MinMQ=20 > filt_complete.vcf
```


Generate statistics for all important items in the vcf files (only works on version 0.1.12b)

```
vcftools --vcf filt_complete.vcf --het --out aaaa
vcftools --vcf filt_complete.vcf --hap-r2 --out aaaa2
vcftools --vcf filt_complete.vcf --freq --out aaaa3
vcftools --vcf filt_complete.vcf --depth --out aaaa4
vcftools --vcf filt_complete.vcf --site-depth --out aaaa5
vcftools --vcf filt_complete.vcf --site-mean-depth --out aaaa6
vcftools --vcf filt_complete.vcf --site-quality --out aaaa7
vcftools --vcf filt_complete.vcf --geno-r2 --out aaaa8
```


Filtering with different MAF

```
vcftools --vcf filt_complete.vcf --maf 0.025 --recode --recode-INFO-all --out filt_completeMAF.025.vcf
vcftools --vcf filt_complete.vcf --maf 0.01 --recode --recode-INFO-all --out filt_completeMAF.01.vcf
vcftools --vcf filt_complete.vcf --maf 0.05 --recode --recode-INFO-all --out filt_completeMAF.05.vcf
```


Used previously to filter 
```
cat filtered_snps.vcf | grep 'PASS\|^#' > filtered_passed_snps.vcf 
```


- MAF? Check literature


**C)    Further filtering + conversion of vcf file. Custom edition**

Now we will locate called SNPs and genotypes and filters on set quality and depth it will
produce an 'N' for every position (SNP/individual) for which no genotype was called by UG
or failing the filter we set.

NOTE: you are now entering in the domain of custom scripts. Although you will be working
with your data more hands-on, you will also produce tables that are custom-made. This
means that there are not necessarily easy format scripts to manipulate the data (as you
need to do this for every-single-program)

```
#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
# minimum genotype quality
my $min_gt_qual = 20;
# minimum mapping quality
my $min_mq = 20;
# minimum SNP quality
my $min_qual = 20;
# minimum depth
my $min_dp = 5; 
# max depth
my $max_dp =240;
            
#old VCF2VERTICAL had format from mpileup: GT:PL:DP:GQ
#THIS VERSION:    has format from GATK-UG: GT:AD:DP:GQ:PL

while(<STDIN>){
        if(eof()){
                print "\n";
        }
        else{
                my $line = "$_";
                chomp $line;
            if($line=~m/^##/){
                    next;
                }
                else{
                        my @fields = split /\t/,$line;
                        my $chrome = shift @fields;
                        my $pos =    shift @fields;
                        my $id =     shift @fields;
                        my $ref =    shift @fields;
                        my $alt =    shift @fields;
                        my $qual =   shift @fields;
                        my $filter = shift @fields;
                        my $info =   shift @fields;
                        my $format = shift @fields;
                        my $mq = 0;
                        if($info=~m/MQ=(\d+)/){
                                $mq = "$1";
                        }
                        my $meta = "$chrome\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format";
                        print "$chrome\t$pos";
                        if($line=~m/^#/){
                                foreach(@fields){
                                        my $long = "$_";
                                        my $name = basename($long,'.bam');
                                        print "\t$name";
                                }
                                print "\n";
                        }
                        else{
                                foreach(@fields){
                                        my $fourbasename = "$_";
                                        my $allele0 = &GT($ref,$alt,$fourbasename);
                                        if(($qual >= $min_qual) && ($mq >= $min_mq)){
                                                print "\t$allele0";


                                        }
                                        else{
                                                print "\tN";
                                        }
                                }
                                print "\n";
                        }
                }
        }
}
sub GT{
        my $ref =shift;
        my $alt =shift;
        my $alt2;
    my $alt3;
    #if there are two alternate alleles:
      if($alt=~m/,/){
                my @alts= split /,/, $alt;
        my $alts_length = @alts;
                $alt=$alts[0];
                $alt2=$alts[1];
        if ($alts_length == 3){
            $alt3=$alts[2];
        }
        }
        my $fourbasename = shift;

        my $substring = "./.";
        if ($fourbasename =~ /\Q$substring\E/){
        return 'N';
    } else {

    my @gtdata = split /:/, $fourbasename;

    #genotype data (PL) is now in the 5th array entry:
        my @genoP = split /,/, $gtdata[4];
    #gq is now 4th:
        my $gq = $gtdata[3];
    #depth is 3rd:
        my $dp = $gtdata[2];
    if ($gq <= $min_gt_qual || $dp <= $min_dp || $dp > $max_dp  ){

                return 'N';
        }
        else{
                my $i =1;
                my $n_match =0;
                my %types = ( 1 => '00', 2 => '01', 3 => '11', 4 => '02', 5 => '12', 6 => '22', 7 => '03', 8 => '13', 9 => '23', 10 => '33');

                my $genotype = 'N';
                #go through each genotype liklihood - they are in order ref/ref, ref/alt1, alt1/alt1, ref/alt2, alt1/alt2, alt2/alt2 
                foreach(@genoP){
                        my $prob = "$_";
                        #PL is L(data given that the true genotype is X/Y) so, bigger is LESS confident
                        if($prob==0){
                                ++$n_match;
                                #Get the genotype based on the position of the 0 
                                if(exists $types{$i}){
                                        $genotype = $types{$i};
                                }
                                else{
                                        $genotype = 'XX';
                                }
                        }
                        ++$i;

                }
        # 00,01,11,02,12,22
        #P(D|CC)=10^{-0.7}, P(D|CA)=1, P(D|AA)=10^{-3.7}, P(D|CG)=10^{-1.3}, P(D|AG)=1e-4 and P(D|GG)=10^{-4.9}.
                #my $genotypequal = $gtdata[2];
                #if($genotypequal<$min_gt_qual){
                #if more than one 0 than return an N
                if($n_match!=1){
                        return 'N';
                }
                else{
                        $genotype =~ s/0/$ref/eg;
                        $genotype =~ s/1/$alt/eg;
                        $genotype =~ s/2/$alt2/eg;
            $genotype =~ s/3/$alt3/eg;
                        $genotype =~ s/\///;
                        if ($genotype eq 'AA'){
                                $genotype = "A";
                        }
                         elsif ($genotype eq 'TT'){
                                $genotype = "T";
                        }
                        elsif ($genotype eq 'CC'){
                                $genotype = "C";
                        }
                        elsif ($genotype eq 'GG'){
                                $genotype = "G";
                        }
                        elsif (($genotype eq 'AC') || ($genotype eq 'CA')){
                                $genotype = "M";
                        }
                        elsif (($genotype eq 'AG') || ($genotype eq 'GA')){
                                $genotype = "R";
                        } 
                        elsif (($genotype eq 'AT') || ($genotype eq 'TA')){
                                $genotype = "W";
                        }
                        elsif (($genotype eq 'CG') || ($genotype eq 'GC')){
                                $genotype = "S";
                        }
                        elsif (($genotype eq 'CT') || ($genotype eq 'TC')){
                                $genotype = "Y";
                        }
                        elsif (($genotype eq 'GT') || ($genotype eq 'TG')){
                                $genotype = "K";
                        }
                        else{
                                $genotype = "N";
                        }
                        return $genotype;
                }
                
                
        }
    }           
                
}          
```

Run as

```
cat filtered_passed_snps.vcf | ~/scripts/vcf2vertical_dep_GATK-UG.pl > snptableUG_pass.tab
```

Get the number of filtered SNPs:

```
cat snptableUG_pass.tab| grep SC2* | wc -l
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

**D)    Filter SNPtable on population genetic parameters**

Get summary statistics and filter SNP table based on minor allele frequency, heterozygosity and missing data using snp_coverage.pl

```
perl /nfs/home/hpcsci/lotteanv/scripts/snp_coverage_p240.pl snptableUG_pass.tab
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
countN<-sapply(lst,FUN=function(x,d){sum(d[,x]=="N",na.rm=T)},d)
```

Calculate the percentage of N per SampleID divided by total number of SNPs

```
percN<-sapply(lst,FUN=function(x,d){sum(d[,x]=="N",na.rm=T)}/nrow(d),d)
```

And make a histogram

```
hist(percN)
```
write.table(percN,"percN_diffMAF.txt",quote=F,sep="\t")
write.table(countN,"countN_diffMAF.txt",quote=F,sep="\t")
```
lst<-gsub(pattern = "X",replacement = "",colnames(d[,-(1:2)]))

Copy these files to local computer --> open new terminal but DO NOT ssh into MCC


Check the number of SNPs gained through each method

```
cat snptableUG.tab.table | wc -l
```



There is another interesting step which can be preformed, which is called phasing. This is where we find out not only what the genotype is at each locus, but also from which of the 2 chromosomes (in a diploid organism) each allele is coming from. This information can be useful when looking at inversions etc

and example of a walkthrough is on https://github.com/mgharvey/seqcap_pop



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