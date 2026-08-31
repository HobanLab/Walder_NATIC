# %%% SCRIPT FOR OBTAINING DENSITY ESTIMATES FOR ALNUS MARITIMA PER STATE ----
library(rFIA)
# Specify the folder in which the CSVs downloaded for each state are being stored
fiaDir <- '/RAID1/FIA_Data/allData' ; setwd(fiaDir)
# Subset entire rFIA database to strictly the states of interest (here, states with known ALMA populations)
ALMA <- readFIA(dir=fiaDir, states=c('GA','OK','MA', 'DE'))
# Clip data to only most recent inventory year, then get density estimate
ALMA_MR <- clipFIA(ALMA, mostRecent = TRUE)
FIAtpa_MR <- tpa(ALMA_MR, totals = TRUE, bySpecies = TRUE)

# Printing out common names, and seaside alder is not included
sort(FIAtpa_MR$COMMON_NAME)
# No alders appear to be included
sort(FIAtpa_MR$SCIENTIFIC_NAME)
