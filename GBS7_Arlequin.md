Arlequin

To transfer structure files to Arlequin, I use PGD Spider
Download, navigate to folder in Terminal and enter:

```
java -Xmx1024m -Xms512m -jar PGDSpider2.jar
```

The conversion script is saved in the PGD Spider dir (Documents/Other). Edit this script for the number of loci. 

Download Arlequin for MacOSx and put the files in the same folder as datafiles

To decide which test to run, open ssdefs.txt and delete the # at the line beginning for every desired test
Note that arl_run.ars is the default setting file. Only edit this when you know what you're doing (I don't)

In terminal, navigate to directory and type

```
./arlecoremac_32bit all_nr.arp
```

OR... I find it easier to run the WinArl version (Get Wine and WineBottler to run it, follow instructions on the website)

Or not...

