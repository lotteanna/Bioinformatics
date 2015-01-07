Reducing size of fragmented reference genome
===


A) In case of fragmented reference, it is better to make a pseudo-scaffold. This is 
a concatenated file with all the contigs togethers, separated with 30 A's.
perl ~/scripts/scaffolds.pl soaprunk61.contig
output files are soaprunk61.contig.pseudo  & scaffold_order.soaprunk61.contig
REMEMBER to translate the SNP table back later to the original genome locations
 with scaffold2contig.v2.pl

B) GATK is not able to handle N's, count and replace these with A's
Count number of N’s in the file as GATK can’t handle N’s
note: whenever counting lines make sure to account for header lines.
Use grep -v '>'  first bascially because you don’t want to count characters in the
header and you should never ever assume the header won't contain some chars.
grep -v is "inverse grep": that is "look for lines without this"
But without accounting for this:

**Option1** (—> gives string of A’s and counts the characters)
tr -cd N < soaprunk61.contig.pseudo | wc -c 
 
*OR*

**Option2**(—> slower option, gives char in lines and then counts the number of lines.
 wc -c won’t work here as there will be twice as many characters, A and newline)
```fgrep -o N soaprunk61.contig.pseudo | wc -l```

In case of N's replace with A's with sed, not done here as wasn't necessary

can also be used to see how many A's were added in the scaffolds.pl script
by counting the difference between soaprunk61.contig and soaprunk61.contig.pseudo
```echo '712102436 - 421159706' | bc```

B) Index .contig (or in case of B, contig.pseudo) file for samtools, bwa and gatk

Check maximum line lengths in file with:
```awk '{ if (x < length()) x = length() } END { print x }' <referencegenome>```
Check for anything shorter that this maximum length (here 70)
```awk 'length($0) < 70' <referencegenome>```
In case not all lines are same length:

Delete white lines
```sed '/^$/d' <filename> > <newfilename>```

-a specifies the indexing algorithm bwa uses. There is another option (IS), for smaller
genomes (<2GB), check bwa man page
```bwa index -a bwtsw <referencegenome>```

```samtools faidx <referencegenome>```