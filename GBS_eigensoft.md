EIGENSOFT
---

Before you do anything! Installing EIGENSOFT on a mac is a pain in the behind!

To transfer STRUCTURE files to EIGENSOFT, I use PGD Spider
Navigate to PDGD Spider folder in Terminal and enter:

```
java -Xmx1024m -Xms512m -jar PGDSpider2.jar
```

Save convert script as STRUCTURE_2_EIGENSOFT in /Documents/Monash/PhD/Analyses/Data
- No phase info
- Individual labels and population identifier present


While converting, got the warnings:
Locus <xxx> has no location information, thus the location is set to chr1.
Not enough memory


Download EIGENSOFT, put it in the desired directory and unpack. Make sure dependencies are present.
The one I found needed to be installed were:
gfortran (https://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
openblas (http://www.openblas.net/)

gfortran is a dmg, so can just be installed by clicking on it
For openblas, navigate to directory and type 

```
make
```

The default directory for installation is /opt/openblas/lib, which is where you want it!
Bluh doesn't work I give up


In EIGENSOFT (terminal based),

```
cd POPGEN
```

../bin/smartpca -p parfile >logfile