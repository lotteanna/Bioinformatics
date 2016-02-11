Preparing reference genome for downstream programs. Also includes a guide on what to do with a highly fragmented reference genome!
===



Load required packages
```
module load samtools
module load bwa
```

Link(s):

http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

---

** Optional) Turn fragmented reference genome into pseudo-scaffold **
In case of fragmented reference, it is better to make a pseudo-scaffold. This is 
a concatenated file with all the contigs togethers, separated with 30 A's.
```
perl ~/scripts/scaffolds.pl soaprunk61.contig
```

output files are soaprunk61.contig.pseudo  & scaffold_order.soaprunk61.contig

REMEMBER to translate the SNP table back later to the original genome locations
 with scaffold2contig.v2.pl

---

**A) Replace N's for A's**

GATK is not able to handle N's, count and replace these with A's
Count number of N’s in the file as GATK can’t handle N’s
note: whenever counting lines make sure to account for header lines.

First, it might be really good practise to understand always what data you're actually looking at. This can be easily done with the option less <filename>

Use grep -v '>' first. Basically because you don’t want to count characters in the
header (starting with the character ">") and you should never ever assume the header won't contain some chars.
grep -v is "inverse grep": that is "look for lines without this"

Without accounting for headers:

**Option1** (—> gives string of A’s and counts the characters)
```
tr -cd N < soaprunk61.contig.pseudo | wc -c 
```

And with accounting for headers:
```
grep -v '>'| tr -cd N < ragweed_17Aug2015_2ABsE.fasta | wc -c 
```

*OR*

**Option2**  —> slower option, gives char in lines and then counts the number of lines.
 wc -c won’t work here as there will be twice as many characters, A and newline

```fgrep -o N soaprunk61.contig.pseudo | wc -l```

In case of N's replace with A's with sed:

can also be used to see how many A's were added in the scaffolds.pl script
by counting the difference between soaprunk61.contig and soaprunk61.contig.pseudo
```echo '712102436 - 421159706' | bc```


In the reference genome there might be use of both capital and non-capital bases. Repeats from RepeatMasker & TandemRepeatFinder are shown in lowercase, non-repeating sequences are shown in upper-case. The term for this is soft-masking. In the opposite hard-masking, these bases will be replaced by N’s. The sequences in lower cases might be lower confidence. Shouldn’t matter as the mapping quality of repetitive sequences is really low, and no SNPs will be called in this region. Options are to leave as is or replace all with the same capital letter (A). Just make sure whatever you decide, that the downstream programs are able to handle lower cases. Just because the downstream programs are always very picky, I'll change all the lowercase with uppercase (might still help in the mapping process to know that a certain part is next to a repetitive region)

the “tr” functions allows for replacements

---

**B) Index** 

Index .contig (or in case of fragmented reference, contig.pseudo) file for samtools, bwa and gatk

Check maximum line lengths in file with:
```
awk '{ if (x < length()) x = length() } END { print x }' <referencegenome>
```

Check for anything shorter that this maximum length (here 70)
```
awk 'length($0) < 70' <referencegenome>
```

In case not all lines are same length:

Delete white lines
```sed '/^$/d' <filename> > <newfilename>```

Otherwise samtools won't be able to index the reference genome and will give the error:
 [fai_build_core] different line length in sequence <Scaffold#>


Program specific index code:

*For bwa*
-a specifies the indexing algorithm bwa uses. There is another option (IS), for smaller
genomes (<2GB), check bwa man page

```
bwa index -a bwtsw <referencegenome>
```

*For samtools*
```
samtools faidx <referencegenome>
```

*For PICARD*
```
java -jar $PICARD/CreateSequenceDictionary.jar REFERENCE=<refgenome> OUTPUT=<refgenome>.dict
```



