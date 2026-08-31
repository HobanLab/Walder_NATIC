# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% EXTRACT DENSITIES FROM WILSON ET AL. 2013 DATASET %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This script reads in, cleans, and processes the data in the Wilson et al. 2013 dataset, and then combines that
# basal area per acre data with data from the FIA dataset to calculate densities for tree species in North America.

# There are 2 sections
# A. DATA CLEANING: This section reads in 3 CSVs from the source data: one for species in the western half
# of the U.S., one for species in the eastern half, and one for species found across the U.S. (at a larger
# resolution of 200km). It creates a new composite dataframe from these 3 CSVs, made up of the nonzero rows
# from each original CSV
# B. DATA ANALYSIS: This section uses the basal area per acre data from the cleaned Wilson et al. 2013 data and
# combines it with the FIA data for the contiguous U.S. to generate density values for different tree species.
pacman::p_load(readr, dplyr, rFIA)
# Specify FIA directory, and directory of Wilson et al. 2013 data
fiaDir <- '/RAID1/FIA_Data/allData' ; setwd(fiaDir)
wilsonDir <- '/home/akoontz/Documents/Indicators/Walder_Indicators/RDS-2013-0013_Data/'

# %%% SECTION A: DATA CLEANNING %%% ----
# Read in the CSVs
western_half <- read_csv(paste0(wilsonDir, "agreement_metrics_western_half.csv"))
eastern_half <- read_csv(paste0(wilsonDir, "agreement_metrics_eastern_half.csv"))
k200 <-  read_csv(paste0(wilsonDir, "agreement_metrics_200k.csv"))
# There are 324 rows in these (note: some of these are "spp.", "Unknown", "NA", etc.)
lapply(list(western_half, eastern_half, k200), nrow)

# A1: CLEAN WESTERN DATA %%%
# Create set of species with nonzero Avg_ba (Average basal area per acre) values in western 
# half and zero Avg_ba values in eastern half
W_species <- western_half %>%
  filter(Avg_ba != 0) %>%
  pull('Common Name')

E_zero_species <- eastern_half %>%
  filter(Avg_ba == 0) %>%
  pull('Common Name')

W_species <- intersect(W_species, E_zero_species)
wilsonWest <- western_half %>% filter('Common Name' %in% W_species)
# Remove these species from eastern half
eastern_half_filtered <- eastern_half %>%
  filter(!'Common Name' %in% W_species)

# A2: CLEAN EASTERN DATA %%%
# Create set of species with nonzero Avg_ba values in eastern half and zero values in western half
E_species <- eastern_half %>%
  filter(Avg_ba != 0) %>%
  pull('Common Name')

W_zero_species <- western_half %>%
  filter(Avg_ba == 0) %>%
  pull('Common Name')

E_species <- intersect(E_species, W_zero_species)
wilsonEast <- eastern_half %>%  filter('Common Name' %in% E_species)
# Remove these species from western half
western_half_filtered <- western_half %>%
  filter(!'Common Name' %in% E_species)

# A3: CLEAN TOTAL DATA %%%
# Create a set of species with nonzero values in 200km dataset, and which aren't present in western/eastern sets
excluded_species <- union(W_species, E_species)
total <- k200 %>% filter(Avg_ba != 0, !'Common Name' %in% excluded_species)

# A4: Combine western, eastern, and total datasets into final dataframe, and remove duplicates
final_df <- bind_rows(wilsonWest, wilsonEast, total)
final_df <- final_df %>%  distinct()
# There are 287 rows in the final dataframe
nrow(final_df)
# Write results to disc
write_csv(final_df, paste0(wilsonDir, 'processedMetrics.csv'))

# %%% SECTION B: DATA ANALYSIS %%% ----
inputData <- paste0(wilsonDir, 'processedMetrics.csv')
usfs <- invisible(read_csv(inputData, col_names = TRUE))
head(usfs)

# Extract FIA data for contiguous U.S.
contigUS <- c("WA","OR","CA","ID","MT","WY","CO","NM","AZ","UT","NV",
              "ND","SD","NE","KS","OK","TX","MN","IA","MO","AR","LA",
              "WI","IL","MI","IN","OH","KY","TN","MS","AL","GA","FL",
              "SC","NC","VA","WV","PA","NY","VT","NH","ME","MA","CT",
              "RI","NJ","DE","MD")
fia <- readFIA(dir = fiaDir, states = contigUS, tables = c("TREE","COND","PLOT"),
               inMemory = FALSE)
# Calculate the Quadratic Mean Diameter (QMD) from the FIA data. This is the diameter of the "average" tree, weighted
# by basal area (meaning larger trees count more than smaller ones)
qmd <- tpa(
  fia,
  bySpecies = TRUE,
  landType = "forest",
  treeType = "live",
  totals = FALSE
) %>%
  mutate(
    QMD = sqrt(BAA / (0.005454 * TPA))
  ) %>%
  select(SPCD, TPA, BAA, QMD)

# The tibble above contains multiple QMD values per species. Before combining these with the basal area per acre values
# provided in the Wilson et al. 2013 dataset (usfs), average QMD values by species
qmd2 <- qmd %>%
  group_by(SPCD) %>%
  summarize(
    TPA = mean(TPA, na.rm=TRUE),
    BAA = mean(BAA, na.rm=TRUE),
    QMD = mean(QMD, na.rm=TRUE)
  )

# Determine how many (and which) species overlap between the FIA dataset and the Wilson et al. 2013 dataset
# Create species lookup table
species_lookup <- distinct(usfs[, c("Spp code", "Scientific Name")])
# Capture overlapping species codes
overlap <- intersect(
  unique(qmd2$SPCD),
  unique(usfs$`Spp code`)
)
# Generate list of overlapping species with names
overlap_species <- species_lookup[species_lookup$`Spp code` %in% overlap,]
overlap_species <- overlap_species[order(overlap_species$`Spp code`),]
nrow(overlap_species)
overlap_species$`Scientific Name`

# Combine USFS and FIA datasets, by matching the species code variables in each
combined <- left_join(usfs,  qmd2,  by = c("Spp code" = "SPCD"))

# First, calculate the basal area per tree by converting the derived QMD value. Then, divide basal area per acre (from Wilson et al. 2013)
# by basal area per tree (from FIA data) to get trees per acre
combined <- combined %>%
  mutate(
    BA_tree = 0.005454 * QMD^2,
    Trees_per_acre = Avg_ba / BA_tree
  )

# Examine outputs: species name, FIA TPA, and TPA calculated based on Wilson et al. 2013 data
finalDensities <- combined[,c(1,2,3,5,41,45)] ; print(finalDensities)

# Write CSV, to output density values per species
write.table(finalDensities, file=paste0(wilsonDir,'calculatedTPAs.csv'),sep=',')
