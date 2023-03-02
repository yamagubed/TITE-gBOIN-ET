################################
###   Software information   ###
################################

All the simulation studies were implemented using R software. The output of sessionInfo() is as follows:

R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8 
[2] LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Iso_0.0-18.1

loaded via a namespace (and not attached):
[1] compiler_4.2.2


####################################################
###   Reproducing the simulation study results   ###
####################################################

1. Run the following programs:

/prog/simulation1/simulation1_gBOIN-ET.R
/prog/simulation1/simulation1_TITE-gBOIN-ET.R
/prog/simulation2/simulation2_gBOIN-ET.R
/prog/simulation2/simulation2_TITE-gBOIN-ET.R

which read (i) /data/true_probbaility.csv for true toxicity and efficacy probabilities, and (ii) /prog/function/myfunction.R containing additional R functions used in each program. The run of the four programs above produces the following csv files containing full simulation studies results:

/data/simulation1_dose_allocation_gBOIN-ET.csv
/data/simulation1_dose_allocation_TITE-gBOIN-ET.csv
/data/simulation1_obd_gBOIN-ET.csv
/data/simulation1_obd_TITE-gBOIN-ET.csv
/data/simulation2_dose_allocation_gBOIN-ET.csv
/data/simulation2_dose_allocation_TITE-gBOIN-ET.csv
/data/simulation2_obd_gBOIN-ET.csv
/data/simulation2_obd_TITE-gBOIN-ET.csv

2. Run the following program:

/prog/tables/tables.R

which produces the following outputs:

/output/Table 1.txt
/output/Table 2.txt
/output/Table S2.txt
/output/Table S3.txt
/output/Table S4.txt
/output/Table S5.txt
/output/Table S6.txt
/output/Table S7.txt

The file names are consistent with the table numbers used in the manuscript.
