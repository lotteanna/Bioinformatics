Creating p-values for pairwise Fst values
===

Code by:
Thibaut Jombart 
Philip Chan 
Lotte van Boheemen 

---

As there seems to be a lack of packages able to calculate p-values for
pairwise FSTs for large datasets, I decided to pull together some code
and get a calculation set up for R. My script is based on a script by Thibaut Jombart at
http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2011-February/000214.html
and all you need it is to use your favourite datafile to read into a genind object
(that is the adegenet package). 

The largest problem I was having was the time it took to do a decent
number of permutations to get the p-values. Doing it on my computer,
it would take ~1 day for 10 perms (one permutation gives one matrix
with pairwise fsts draws from a random distribution, your 0
hypothesis), but you need 1000 (=1000 pfst matrices) at least (I do
have 85 pops so this is contributing to the duration). For this
reason, what you can do is send it to the cluster in an array, where
each array is analysing a subset (let's say 10 matrices). Then, all
you need to do is concatenate the files and do an additional step
which actually calculates your p-values. What that final step does is
testing if your observed fst matrix is within the 1000 random
matrixes, for every value within the matrix. If an observed value (say
pop 1 vs pop 2) is greater then 950 of the 1000 values of pop1 vs pop2
in the permutated matrices, it is significant (p <0.05).

**This script DOES NOT work on the MCC, you need to load your data onto Monarch.**

```
##R-script to draw perutated values from the pairwise matrix on the cluster

library("adegenet")
library("hierfstat")
library("ggplot2")

#this section is my way of transforming a structure file into genind
allr<-read.structure("allcra473_240_in.stru", n.ind=473, n.loc=1000, onerowperind=FALSE, col.lab=1, col.pop=2, row.marknames=1,NA.char="-9", pop=85, ask=FALSE, quiet=FALSE)

#do the calculation
mat.perm <- pairwise.fst(allr, pop=sample(pop(allr)), res.type="matrix")
#write the output to a file with taskID as identifier
write.table(mat.perm,paste("perm_allr.txt"),sep="\t")
```


Below is the job script. in this script, just change the time  and in one of your final lines change

```cp -p ../../all* . ``` to ```
cp -p ../../<whatever your structure file is named> . ```


```
#!/bin/env bash
#SBATCH --job-name=matperm
#SBATCH --time=50:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=matperm.out
#SBATCH --array=1-1000

module load R/3.2.2
ID=${SLURM_ARRAY_TASK_ID}

P=`expr $ID / 100`
C=`expr $ID % 100`

PDIR=`printf "%03d" $P`
CDIR=`printf "%03d" $C`

echo "Working on $PDIR / $CDIR"
if [ ! -d $PDIR ]
then
  mkdir $PDIR
fi

cd $PDIR
if [ ! -d $CDIR ]
then
  mkdir $CDIR
fi

cd $CDIR
IND=`expr $ID - 1`
INDEX=`expr $IND \* 2`

if [ -f bf_$INDEX.bf ]
then
  echo "Output file for $PDIR / $CDIR already present!"
  exit 0
fi

pwd

cp -p ../../all* .
cp -p ../../matperm.R .

ulimit -s 100000

Rscript matperm.R 
```


Next job is to past all the matrices together. Below code will search
within your parent dirs and child dirs within and concatenates
everything starting with perm (so all the permutated matrices)
together into 1 file. This can be done in a job, so I haven't put all
the R code in one block, as for me it runs too long within the
console. Replace the path name with your dir (you get this in R with
```getwd()```), and replace perm_allr.txt with whatever name you gave
each permutated matrix:


```
module load R
R
dat.files  <- list.files(path="/mnt/lustre/projects/p2016050002/pfstat/",
                                recursive=T,
                                pattern="perm_allr.txt"
                                ,full.names=T)
readDatFile <- function(f) {
  dat.fl <- read.table(f) 
}
mat.perm<- lapply(dat.files, readDatFile)
write.table(mat.perm,"mat.perm")
```

Then, you have to get the actual pairwise FST values for your dataset (can be done in a job again):

```
module load R
R
mat.obs <- pairwise.fst(allr, res.type="matrix") 
write.table(mat.obs,"mat.obs")
```

Finally, you can compare between your actual matrix and all the permutated matrices:

```
library(ade4)

mat.obs<-read.table("mat.obs")
mat.perm<-read.table("mat.perm")

allTests <- list()
 for(i in 1:(nrow(mat.obs)-1)){
   for(j in 2:nrow(mat.obs)){
   allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:1000, 
   function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")
   }
}
```

In the console, save the RData for later use:
```
dput(allTests,"allTests")
```

Download to the local machine and open in a R Markdown script

```
load("allTests")
```

You can get the basic information as follows

```
ls.str((allTests)
```

Knit the Markdown script and copy-paste to a text editor. In terminal, delete all the
non-sense information:

```
egrep -v 'sim|obs|alter|rep|expvar|cal' allTests.txt > redTests.txt
```

Remove all the other unnecessary stuff (in the most clunky way)
```
sed -i -e 's/## //' redTests.txt
```

