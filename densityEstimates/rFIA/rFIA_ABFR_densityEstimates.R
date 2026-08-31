# %%% SCRIPT FOR OBTAINING DENSITY ESTIMATES FOR ABIES FRASERI PER STATE ----
library(rFIA)
# Specify the folder in which the CSVs downloaded for each state are being stored
fiaDir <- '/RAID1/FIA_Data/allData' ; setwd(fiaDir)
# Subset entire rFIA database to strictly the states of interest (here, states with known ABFR populations)
ABFR <- readFIA(dir=fiaDir, states=c('GA','TN','NC', 'VA', 'WV'))
# Clip data to only most recent inventory year, then get density estimate
ABFR_MR <- clipFIA(ABFR, mostRecent = TRUE)
# FIAtpa_MR <- tpa(ABFR_MR, totals = TRUE, bySpecies = TRUE)
# FIAtpa_MR[which(ABFR_MR$COMMON_NAME=='Fraser fir'),5]

# # Georgia density
# ABFR_GA <- readFIA(dir=fiaDir, states=c('GA'))
# # Clip data to only most recent inventory year, then get density estimate
# ABFR_GA <- clipFIA(ABFR_GA, mostRecent = TRUE)
# ABFR_GA_TPA <- tpa(ABFR_GA, totals = TRUE, bySpecies = TRUE)
# GA_TPA <- as.numeric(ABFR_GA_TPA[which(ABFR_GA_TPA$COMMON_NAME=='Fraser fir'),5])

# Tennessee density
ABFR_TN <- readFIA(dir=fiaDir, states=c('TN'))
# Clip data to only most recent inventory year, then get density estimate
ABFR_TN <- clipFIA(ABFR_TN, mostRecent = TRUE)
ABFR_TN_TPA <- tpa(ABFR_TN, totals = TRUE, bySpecies = TRUE)
TN_TPA <- as.numeric(ABFR_TN_TPA[which(ABFR_TN_TPA$COMMON_NAME=='Fraser fir'),5])
# Calculate density in trees per km2, and report both trees per acre and trees per km2
TN_TPKM <- TN_TPA/0.00404686 ; paste0('TN TPA: ', TN_TPA, '; TN TPKM: ', TN_TPKM)

# North Carolina density
ABFR_NC <- readFIA(dir=fiaDir, states=c('NC'))
# Clip data to only most recent inventory year, then get density estimate
ABFR_NC <- clipFIA(ABFR_NC, mostRecent = TRUE)
ABFR_NC_TPA <- tpa(ABFR_NC, totals = TRUE, bySpecies = TRUE)
NC_TPA <- as.numeric(ABFR_NC_TPA[which(ABFR_NC_TPA$COMMON_NAME=='Fraser fir'),5])
# Calculate density in trees per km2, and report both trees per acre and trees per km2
NC_TPKM <- NC_TPA/0.00404686 ; paste0('NC TPA: ', NC_TPA, '; NC TPKM: ', NC_TPKM)

# Virginia density
ABFR_VA <- readFIA(dir=fiaDir, states=c('VA'))
# Clip data to only most recent inventory year, then get density estimate
ABFR_VA <- clipFIA(ABFR_VA, mostRecent = TRUE)
ABFR_VA_TPA <- tpa(ABFR_VA, totals = TRUE, bySpecies = TRUE)
VA_TPA <- as.numeric(ABFR_VA_TPA[which(ABFR_VA_TPA$COMMON_NAME=='Fraser fir'),5])
# Calculate density in trees per km2, and report both trees per acre and trees per km2
VA_TPKM <- VA_TPA/0.00404686 ; paste0('VA TPA: ', VA_TPA, '; VA TPKM: ', VA_TPKM)

# Average of TN and NC TPKMs
mean(TN_TPKM, NC_TPKM)

# Sum of population sizes
12368+12538+22149+15198+11329+160

# # West Virginia density
# ABFR_WV <- readFIA(dir=fiaDir, states=c('WV'))
# # Clip data to only most recent inventory year, then get density estimate
# ABFR_WV <- clipFIA(ABFR_WV, mostRecent = TRUE)
# ABFR_WV_TPA <- tpa(ABFR_WV, totals = TRUE, bySpecies = TRUE)
# WV_TPA <- as.numeric(ABFR_WV_TPA[which(ABFR_WV_TPA$COMMON_NAME=='Fraser fir'),5])
