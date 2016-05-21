
And below is the R script, called matperm.R in the jobscript in the previous email. You need to check the "read.structure" part of the code to be correct for your infile, i.e. define n.ind (number of individuals in your infile), number of lines per individual (if one, onerowperind=TRUE), specify column with individual names (col.lab=<column number>), names of populations (col.pop=<column number>), NA character used and number of populations in infile (pop=<#>). Best is to read it into R on your computer before submitting, and check if it understands what your are feeding it. If it tells you below, you are doing great
 Converting data from a STRUCTURE .stru file to a genind object... 

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

