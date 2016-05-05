**TESS**
===

http://membres-timc.imag.fr/Olivier.Francois/tess.html

Written by Alexandra Pavlova (2011) & Lotte van Boheemen (2016)


**A) Create TESS files**

Exact same as STRUCTURE, except it has to contain longitude and latitude (in that order) *before* the genotype data. This can be either 1 row per individual (so 2 columns per locus, as we're dealing with a diploid organism) or 2 rows per individual. As we are following the format of the structure file as produced in GBS6_structure, below code will assume data is organised as 2 rows per individual (default setting, otherwise use -spy)

Missing data is noted as -9. Preferably coordinates need to be different per individual, but as far as I can tell, this is not necessary.

Example input file

```
SampleID	pop_no	long	lat	SC2ABSE_564_2_21120	SC2ABSE_398_3_100654	SC2ABSE_398_3_100779
1-010908-1-12	1	18.7826	47.642	-9	-9	4
1-010908-1-12	1	18.7826	47.642	-9	-9	4
1-010908-1-27	1	18.7826	47.642	-9	-9	-9
1-010908-1-27	1	18.7826	47.642	-9	-9	-9
1-160808-1-1	2	8.94937	46.0325	-9	-9	3
1-160808-1-1	2	8.94937	46.0325	-9	-9	4
```

Extra columns in the datafile (e.g. sampleID and pop_no) are indicated by -c (so for 2 additional columns, this is -c2).
Additional rows (e.g. locus names) are indicated by -r (so for 1 additional row, this is -r1)

Be careful when editing these files in Excel, it might not be saved in the correct format (e.g. all in one line)


**B) Prepare TESS run for cluster**

This is specifically designed to set up a job and submit it to the Monash Campus Cluster (MCC)


First, make the script start.sh

```
#start.sh
#!/bin/sh -l
#$ -S /bin/sh
#$ -cwd

RTYPE=0
FILENAME="!"
KPARM=0
MAX_K=0
MAX_QSUB=0

while [ ! -e $FILENAME ]; do
  echo -e "Input filename: \c"
  read FILENAME
done  
echo -e "Minimum K parameter: \c"
read KPARM
echo -e "Maximum K parameter: \c"
read MAX_K
echo -e "Number of replicates: \c"
read MAX_QSUB
echo -e "CAR or BYM"
read RTYPE

if [ $RTYPE = "CAR" ]; then
  while [ $KPARM -le $MAX_K ]; do
    CMD=`qsub -t 1-$MAX_QSUB arrayCAR.job $KPARM $FILENAME`
    echo "Submitted k=$KPARM"
    let KPARM=KPARM+1
  done
else
  if [ $RTYPE = "BYM" ]; then  
    while [ $KPARM -le $MAX_K ]; do
      CMD=`qsub -t 1-$MAX_QSUB arrayBYM.job $KPARM $FILENAME`
      echo "Submitted k=$KPARM"
      let KPARM=KPARM+1
    done
  fi    
fi
echo "Done"
```




====

Below scripts will create directories for each run and give a walkthrough how to download and subsequently use TESS gui to manipulate and look at the data.

```
#written by Alexandra Pavlova 2011
#arrayCAR.job
#!/bin/sh -l
#$ -S /bin/sh
#$ -cwd
#$ -l h_rt=6:00:00
#$ -l h_vmem=4G

#First parameter - K parameter
#Second parameter - Input file name

#Calculate directory name for the run
DIRNAME=$2_CAR_K$1_RUN_$SGE_TASK_ID

#Create directory
mkdir $DIRNAME

#copy input file
cp $2 $DIRNAME/$2

#Load Tess
module load tess/2.3


#CHANGE THE FOLLOWING PARAMETERS:
#-N: NUMBER OF INDIVIDUALS
#-A: Ploidy (1 = Haploid, 2 = Diploid, ...)
#-L: NUMBER OF LOCI
#-P: interaction parameter of HMRF
#-T: Degree of Trend
#-D: Dirichlet allele frequency model
#-r: NUMBER OF ADDITIONAL ROWS (E.G. LOCUS NAMES)
#-c: NUMBER OF ADDITIONAL COLUMNS (E.G. SAMPLE AND POPULATION NAMES
#-S: NUMBER OF ITERATIONS
#-B: BURNIN
#-upy: update psi (upy = yes, upn = no (default)

#Execute
CMD=`tessmV11CAR -F$2 -N41 -A2 -L14 -K$1 -P0.6 -T2 -D1.0 -r1 -c2 -S100000 -B30000 -upy -i$DIRNAME -o$DIRNAME`
```



arrayBYM.job [Edit parameters after #Execute as per manual EXCEPT –F, –K, -i and –o: they must stay –F$2, –K$1, -i$DIRNAME,  -o$DIRNAME. You must change –N (number of individuals) and -L (number of loci)]

```
#arrayBYM.job
#written by Alexandra Pavlova 2011
#!/bin/sh -l
#$ -S /bin/sh
#$ -cwd

#$ -l h_rt=6:00:00
#$ -l h_vmem=4G

#First parameter - K parameter
#Second parameter - Input file name

#Calculate directory name for the run
DIRNAME=$2_BYM_K$1_RUN_$SGE_TASK_ID

#Create directory
mkdir $DIRNAME

#copy input file
cp $2 $DIRNAME/$2

#Load Tess
module load tess/2.3

#Execute
CMD=`tessmV11BYM -F$2 -N251 -A2 -L14 -K$1 -T2 -D1.0 -r1 -c2 -spy -S100000 -B30000 -i$DIRNAME -o$DIRNAME`
```

---


**C)    Submit to cluster**

Call start.sh

```
bash start.sh
```

When asked input filename enter filename: <filename>
When asked min K enter minimum number of K (usually 2): <2>
When asked max K enter maximum number of clusters you think you have: <10> or whatever
When asked Number of replicates: <100> 100 was recommended by authors
When asked CAR or BYM: <CAR> or <BYM>, these are two admixture models in TESS


After this, your jobs are submitted and result folders will be put in the same directory where your input file was. Check their status with <qstat> 

---

**INTERMEZZO**

If you're interested, like I was, I will explain what all these things actually do

Starting with start.sh (ooo it's Friday afternoon and I cannot get punnier).

The core of this script is:

```
qsub -t 1-$MAX_QSUB arrayCAR.job $KPARM $FILENAME`
```

What this comes down to is that with a max number of reps of 10, and a K form 2-4 this corresponds with:

```
qsub -t 1-10 arrayCAR.job 2 au41
qsub -t 1-10 arrayCAR.job 3 au41
qsub -t 1-10 arrayCAR.job 4 au41
```

In arrayCAR.job it then calls

```
DIRNAME=$2_CAR_K$1_RUN_$SGE_TASK_ID
```

In above example this means:

```
DIRNAME=au41_CAR_K2_1_id
```

So below means what?

```
mkdir $DIRNAME
cp $2 $DIRNAME/$2
CMD=`tessmV11CAR -F$2 -N41 -A2 -L1000 -K$1 -P0.6 -T2 -D1.0 -r1 -c2 -S1000 -B300 -upy -i$DIRNAME -o$DIRNAME`
```

That's right, we are making a dir carrying the name of every run for every K, and tessCAR will read and write from and to this newly created dir.

---

**D)    Post-run**

Copy all directories with results to your computer (it’s going to be (#replicatesX(max K-1)) of
them, so 900 in the example above).


2. You could fool TESS GUI into believing that these results were generated via GUI. a. Create a
project in TESS (preferably with the same name as your input file). Do NOT use “_CAR_” ,” _BYM_” or
“_K” as any part of the file name for files in that target directory (except for folders created by
Tess) as these are used by apvtessext.exe to read model and Kmax into your new project file. b. Move
all your results folders into the same directory where your input file is c. Run a little program
apvtessext.exe (written by Anton Polesskiy) to update TESS GUI project file (.tp) with information
about your runs. This program will ask you to change some running parameters. You could leave them
as is at your own risk. If you leave the “admixture mode” on “auto”, the program will update the
project file (by reading CAR or BYM from folder names into your project file). In addition, two
project files will be created- one for all CAR and one for all BYM runs.If you only have CAR or BYM
you may choose these options (or use “auto”). d. Run TESS GUI as you would normally do to analyse
your results (follow the TESS manual for guidance). TESS could create the file for CLUMPP for you,
so you could average the results from, say, 10% of your best runs for a given K. Also, you could map
your results of admixture on the nice looking map using R script (explained in TESS manual).

===


arrayCAR.job [Edit parameters after #Execute as per manual EXCEPT –F, –K, -i and –o: they must stay
–F$2, –K$1, -i$DIRNAME,  -o$DIRNAME. You must change –N (number of individuals) and -L (number of
loci)]

```
#arrayCAR.job

#!/bin/sh -l
#$ -S /bin/sh
#$ -cwd

#$ -l h_rt=6:00:00
#$ -l h_vmem=4G

#First parameter - K parameter
#Second parameter - Input file name

#Calculate directory name for the run
OUTPUTNAME=$2_CAR_K$1_RUN_$SGE_TASK_ID

#copy input file
cp $2 $OUTPUTNAME

#Load Tess
module load tess/2.3


#CHANGE N and L!
#-N: NUMBER OF INDIVIDUALS
#-A: Ploidy (1 = Haploid, 2 = Diploid, ...)
#-L: NUMBER OF LOCI
#-P: interaction parameter of HMRF
#-T: Degree of Trend
#-D: Dirichlet allele frequency model
#-r: NUMBER OF ADDITIONAL ROWS (E.G. LOCUS NAMES)
#-c: NUMBER OF ADDITIONAL COLUMNS (E.G. SAMPLE AND POPULATION NAMES
#-S: NUMBER OF ITERATIONS
#-B: BURNIN
#-upy: update psi (upy = yes, upn = no (default)


#Execute
CMD=`tessmV11CAR -F$OUTPUTNAME -N41 -A2 -L1000 -K$1 -P0.6 -T2 -D1.0 -r1 -c2 -S1000 -B300 -upy`
```

---


**INTERMEZZO**

If you're interested, like I was, I will explain what all these things actually do

Starting with start.sh (ooo it's Friday afternoon and I cannot get punnier).

The core of this script is:

```
qsub -t 1-$MAX_QSUB arrayCAR.job $KPARM $FILENAME`
```

What this comes down to is that with a max number of reps of 10, and a K form 2-4 this corresponds with:

```
qsub -t 1-10 arrayCAR.job 2 au41
qsub -t 1-10 arrayCAR.job 3 au41
qsub -t 1-10 arrayCAR.job 4 au41
```

In arrayCAR.job it then calls

```
OUTPUTNAME=$2_CAR_K$1_RUN_$SGE_TASK_ID
```

In above example this means:

```
OUTPUTNAME=au41_CAR_K2_runnumber
```

So below means what?

```
cp $2 $OUTPUTNAME
CMD=`tessmV11CAR -F$OUTPUTNAME -N41 -A2 -L1000 -K$1 -P0.6 -T2 -D1.0 -r1 -c2 -S1000 -B300 -upy `
```

That's right, we are making a new filename carrying the name of every run for every K, and copying all the data into this. TessCAR will read and write from and to this newly created filename.

---

**D)    Post-run**


You are on your own with your results, as per TESS manual.

