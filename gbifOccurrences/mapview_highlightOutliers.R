# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% MAPVIEW: HIGHLIGHT OUTLIER OCCURRENCES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This script maps coordinates using mapview and highlights potential outliers
#     1. Spatial outliers
#     2. Temporal outliers
#     3. Attributes: outlier institution, basisOfRecord, taxonRank, coordinateUncertainty
#     4. Points in urbanized areas
# The tukeyFlag function flags continuous variables (spatial and temporal outliers; 
# the FenceMultiplier variables set the threshold of how unique values need to be
# to be flagged as outliers.
# The rareFlag function flags categorical variables (attributes), and the rareProp
# variable sets the proportion that determines whether an occurrence is considered
# an outlier or not.
library(sf)
library(mapview)
library(rnaturalearth)  # For urban area polygons (Natural Earth, 1:10m)

# READ IN CSV ----
# Specify the file path to the relevant CSV below
csvFile <- 
  '/home/akoontz/Documents/Indicators/Walder_Indicators/Scripts/GBIF_occurrences/2026-08-19_NATIClist/species_csvs/Austin/Torreya_californica_1783n_2026-08-19.csv'
pts <- read.csv(file=csvFile, header=TRUE)  
names(pts)  # Check names of columns with the lat/long

# CONVERT TO SF OBJECT ----
# Need to specify a coordinate reference system (crs; WGS84 = 4326; NAD83 = 4269)
pts_sf <- st_as_sf(pts, coords = c("decimallongitude", "decimallatitude"), crs = 4326)

# OUTLIER SETTINGS AND HELPER FUNCTIONS ----
# Tukey fence: values beyond Q1 - multiplier*IQR or Q3 + multiplier*IQR are outliers
# larger multipliers = stricter fence = fewer records flagged
spatialFenceMultiplier <- 3.5
temporalFenceMultiplier <- 1.5
# categorical values making up less than this proportion of records are outliers
rareProp <- 0.01

# Function for flagging continuous variables (spatial/temporal outliers)
tukeyFlag <- function(x, fenceMultiplier = 1.5, lower = TRUE, upper = TRUE){
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  out <- rep(FALSE, length(x))
  if(upper) out <- out | (x > q[2] + fenceMultiplier*iqr)
  if(lower) out <- out | (x < q[1] - fenceMultiplier*iqr)
  out[is.na(out)] <- FALSE
  return(out)
}

# Function for flagging categorical variables attributes)
rareFlag <- function(x, prop = rareProp){
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "missing"
  freq <- table(x) / length(x)
  return(x %in% names(freq)[freq < prop])
}

# 1. SPATIAL OUTLIERS ----
# Nearest-neighbor distance (km) calculated among unique locations,
# so duplicate records at the same coordinates don't mask isolation
coords <- st_coordinates(pts_sf)
locKey <- paste(coords[,1], coords[,2])
uniqLoc <- pts_sf[!duplicated(locKey), ]
distMat <- units::drop_units(st_distance(uniqLoc))
diag(distMat) <- Inf
nnDist <- apply(distMat, 1, min) / 1000
pts_sf$nnDist_km <- nnDist[match(locKey, locKey[!duplicated(locKey)])]
rm(distMat)
# Flag unusually large nearest-neighbor distances (upper Tukey fence)
pts_sf$flag_spatial <- tukeyFlag(pts_sf$nnDist_km, fenceMultiplier = spatialFenceMultiplier, lower = FALSE)
# Alternative: flag the uppermost quartile of nearest-neighbor distances
# pts_sf$flag_spatial <- pts_sf$nnDist_km > quantile(pts_sf$nnDist_km, 0.75)

# 2. TEMPORAL OUTLIERS ----
pts_sf$flag_temporal <- tukeyFlag(pts_sf$year, fenceMultiplier = temporalFenceMultiplier)

# 3. ATTRIBUTE OUTLIERS ----
# Rare institution, basisOfRecord, or taxonRank values
pts_sf$flag_institution <- rareFlag(pts_sf$institutionCode)
pts_sf$flag_basisOfRecord <- rareFlag(pts_sf$basisOfRecord)
pts_sf$flag_taxonRank <- rareFlag(pts_sf$taxonRank)
# Unusually high coordinate uncertainty (log scale, since values are highly skewed)
pts_sf$coordUnc_m <- suppressWarnings(as.numeric(pts_sf$coordinateUncertaintyInMeters))
pts_sf$flag_coordUnc <- tukeyFlag(log10(pts_sf$coordUnc_m + 1), lower = FALSE)
pts_sf$flag_attribute <- pts_sf$flag_institution | pts_sf$flag_basisOfRecord |
  pts_sf$flag_taxonRank | pts_sf$flag_coordUnc

# 4. OCCURRENCES IN URBAN AREAS ----
sf_use_s2(FALSE)  # Avoids errors from invalid geometries in the urban polygons
urban <- ne_download(scale = 10, type = "urban_areas", category = "cultural", returnclass = "sf")
urban <- st_make_valid(st_transform(urban, 4326))
urban <- st_crop(urban, st_bbox(pts_sf))
pts_sf$flag_urban <- lengths(st_intersects(pts_sf, urban)) > 0
sf_use_s2(TRUE)

# SUMMARIZE FLAGS ----
flagCols <- c("flag_spatial", "flag_temporal", "flag_institution", "flag_basisOfRecord",
              "flag_taxonRank", "flag_coordUnc", "flag_urban")
flagDF <- st_drop_geometry(pts_sf)[, flagCols]
pts_sf$flag_any <- rowSums(flagDF) > 0
pts_sf$flagReason <- apply(flagDF, 1, function(r) paste(sub("flag_", "", flagCols)[r], collapse = "; "))
colSums(st_drop_geometry(pts_sf)[, c(flagCols, "flag_any")])  # number of records per flag

# CREATE MAP ----
# Each flag type is its own layer, toggled on/off with the widget on the left
flagLayer <- function(col, name, color){
  sub <- pts_sf[pts_sf[[col]], ]
  if(nrow(sub) == 0) return(NULL)
  mapview(sub, layer.name = name, col.regions = color, color = color, alpha.regions = 0.8, cex = 6)
}
# Create basemap, then add layers
baseMap <- mapview(pts_sf, layer.name = "All occurrences", col.regions = "grey70", color = "grey30",
                   alpha.regions = 0.5, cex = 4, map.types = c("Esri.WorldImagery", "Esri.WorldTopoMap"))
layers <- list(
  if(nrow(urban) > 0) mapview(urban, layer.name = "Urban areas", col.regions = "yellow", alpha.regions = 0.2),
  flagLayer("flag_spatial",   "Spatial outliers",   "red"),
  flagLayer("flag_temporal",  "Temporal outliers",  "orange"),
  flagLayer("flag_institution",   "Rare institution",       "purple"),
  flagLayer("flag_basisOfRecord", "Rare basisOfRecord",     "magenta"),
  flagLayer("flag_taxonRank",     "Rare taxonRank",         "green"),
  flagLayer("flag_coordUnc",      "High coord uncertainty", "blue"),
  flagLayer("flag_urban",     "Urban occurrences",  "cyan")
)
outlierMap <- Reduce(`+`, Filter(Negate(is.null), layers), baseMap)
outlierMap
# Click on each point to see its attributes, including nnDist_km and flagReason

# SAVE OUT (STEPS WHEN USING RSTUDIO) ----
# in figure panel, choose 'export' -> 'save as web page' ; this creates an html file
# with the html file downloaded locally (i.e., does not work if file is in an online shared folder), 
# open the html file in a browser. This is your interactive map!
