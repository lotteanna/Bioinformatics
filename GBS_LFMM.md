LFMM

**)Impute missing data

LFMM is able to handle missing data itself ( set -m to TRUE), but it's adviced to impute >10% missing genotypes with IMPUTE2 or MENDEL-IMPUTE

**) prepare environmental data**

It is adviced to run LFMM on a summary of the environmental variables, e.g. retrieved with PCA. Can also run the variables all together. If -a is set to TRUE, LFMM will run all environmental variables simultaneously

**) format datafile**

Can be formatted from PED (which can be retrieved from PGDspider) and is build-in in the LFMM program:

```
ped2lfmm examples/format_example/example.ped
```

**)Show results**

To display a manhattan plot for the results, consider the R package qqman.