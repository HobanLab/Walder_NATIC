# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GBIF occurrence download + CoordinateCleaner QC pipeline %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This script downloads occurrences from GBIF (using rgbif) and then
# cleans those occurrences (using the CoordinateCleaner package). Every 
# stage prints a timestamped message (visible in run_logtxt), and the
# two slow steps (GBIF download prep and cleaning using Coordinate Cleaner) are
# checkpointed to disk as .rds files. If the script dies partway through and
# you re-run it, it will skip any stage whose checkpoint file already exists,
# rather than re-downloading or re-cleaning from scratch. Delete a checkpoint
# file (or set FORCE_RERUN <- TRUE below) to force that stage to redo.
#
# Steps:
#   1. Read a species list, match to GBIF backbone accepted taxon keys.
#   2. Enumerate synonyms per species (documentation only -- GBIF's taxonKey
#      filter already includes synonym-filed occurrences automatically).
#   3. Submit one occ_download() covering all species, restricted to:
#        - North America (US, Canada, Mexico, Caribbean)
#        - years 1975-2025
#        - hasCoordinate = TRUE, hasGeospatialIssue = FALSE
#   4. Run the CoordinateCleaner record-level pipeline (cc_val ... cc_dupl),
#      a coordinate-precision filter, and the two dataset-level tests
#      (cd_ddmm, cd_round). cc_outl, cd_ddmm, and cd_round are parallelized
#      (see PARALLEL PLAN below)
#
# Requires a free GBIF account (https://www.gbif.org/user/profile) for
# occ_download(). Store credentials once in your .Renviron (recommended):
#
#   usethis::edit_r_environ()
#   # then add these three lines, save, and restart R:
#   GBIF_USER=your_username
#   GBIF_PWD=your_password
#   GBIF_EMAIL=your_email@example.com
# ==============================================================================
pacman::p_load(rgbif, taxize, CoordinateCleaner, dplyr, countrycode, 
               readr, purrr, tibble, future, furrr)

# %%% 1. Parameters ---------------------------------------------------------------
setwd('/home/akoontz/Documents/Indicators/Walder_Indicators/Scripts/GBIF_occurrences/')
species_file <- "speciesLists/SpeciesList_2026-08-20_40n_MexicanSpecies.csv"   # one species name per row, no header
out_dir      <- "gbif_clean_output"
dir.create(out_dir, showWarnings = FALSE)

year_range <- c(1975, 2025)

# North America: US, Canada, Mexico + Caribbean (ISO 3166-1 alpha-2 codes)
countries <- c(
  "US", "CA", "MX",                                            # US, Canada, Mexico
  "AG", "BS", "BB", "CU", "DM", "DO", "GD", "HT", "JM", "KN",   # Caribbean nations
  "LC", "VC", "TT", "AI", "AW", "BQ", "VG", "KY", "CW", "GP",   # + territories
  "MQ", "MS", "PR", "SX", "TC", "VI"
)
# Max acceptable coordinate uncertainty, in km. For range-scale/species-level work,
# a much tighter value like 1-10 km is usually more appropriate.
coord_uncertainty_km <- 10

# ---- PARALLEL PLAN ---------------------------------------------------------------
# Specify number of workers (low number, to leave space on the server)
N_WORKERS <- min(18, future::availableCores() - 2)
plan(multicore, workers = N_WORKERS)   # fork-based parallelization; Linux only
FORCE_RERUN <- FALSE   # set TRUE to ignore existing checkpoints and redo everything

# ---- Logging helper ---------------------------------------------------------------
log_msg <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}
script_start <- Sys.time()
pdf(NULL)   # suppress default Rplots.pdf when running non-interactively
log_msg("=== Pipeline started (", N_WORKERS, " parallel workers) ===")

# %%% 2. Read species list -----------------------------------------------------
species_list <- read_csv(species_file, col_names = "species", show_col_types = FALSE) %>%
  pull(species) %>%
  trimws() %>%
  unique()
log_msg(length(species_list), " species read from ", species_file)

# %%% 3. Match names to GBIF backbone accepted taxon keys ----------------------
# GBIF's own `taxonKey` occurrence filter already includes synonym-filed records 
# automatically -- per GBIF's API documentation, "All included and synonym taxa 
# are included in the search" for a taxonKey query. So once a species is matched 
# to its accepted GBIF backbone usageKey, the occ_download() below will already 
# capture records filed under older/alternate names without any extra step. 
# Synonyms are still enumerated explicitly below (via GBIF's own backbone, 
# not an external service) purely to have a documented record of what was folded in.
log_msg("Matching ", length(species_list), " names to GBIF backbone...")

backbone <- name_backbone_checklist(species_list)

# If there are any unmatched names, write these to a CSV
unmatched <- backbone %>% filter(matchType %in% c("NONE", "HIGHERRANK"))
if (nrow(unmatched) > 0) {
  log_msg(nrow(unmatched), " names did not cleanly match the GBIF backbone; see unmatched_names.csv")
  write_csv(unmatched, file.path(out_dir, "unmatched_names.csv"))
}

accepted <- backbone %>% filter(matchType %in% c("EXACT", "FUZZY"))

log_msg(nrow(accepted), " of ", length(species_list),
        " species matched to a GBIF backbone accepted usage")

# %%% 4. Enumerate synonyms for documentation ------------------------------------
# Primary source: GBIF's own backbone. Optionally cross-checked against Kew's 
# Plants of the World Online via taxize, which is wrapped in tryCatch 
# since that external API is known to rate-limit/block programmatic access 
# (HTTP 403/429/502). Parallelized.
get_gbif_synonyms <- function(key) {
  tryCatch({
    syn <- name_usage(key = key, data = "synonyms")$data
    if (is.null(syn) || nrow(syn) == 0) return(character(0))
    unique(syn$canonicalName[!is.na(syn$canonicalName)])
  }, error = function(e) character(0))
}

get_pow_synonyms_optional <- function(sp) {
  tryCatch({
    syn_df <- taxize::synonyms(sp, db = "pow")[[1]]
    if (is.null(syn_df) || !is.data.frame(syn_df) || nrow(syn_df) == 0) return(character(0))
    name_col <- intersect(c("name", "fullName", "scientificName"), names(syn_df))[1]
    if (is.na(name_col)) return(character(0))
    unique(as.character(syn_df[[name_col]]))
  }, error = function(e) character(0))
}

log_msg("Enumerating synonyms for ", nrow(accepted), " species (", N_WORKERS, " workers)...")

synonym_table <- accepted %>%
  select(species = verbatim_name, usageKey) %>%
  mutate(
    gbif_synonyms = future_map(usageKey, get_gbif_synonyms, .progress = FALSE),
    pow_synonyms  = future_map(species, get_pow_synonyms_optional, .progress = FALSE)
  )

synonym_report <- synonym_table %>%
  mutate(
    gbif_synonyms = map_chr(gbif_synonyms, ~ paste(.x, collapse = "; ")),
    pow_synonyms  = map_chr(pow_synonyms, ~ paste(.x, collapse = "; "))
  )
write_csv(synonym_report, file.path(out_dir, "synonym_mapping.csv"))

gbif_keys <- unique(accepted$usageKey)
log_msg(length(gbif_keys), " unique GBIF accepted taxon keys will be queried ",
        "(synonym occurrences are included automatically by GBIF's taxonKey filter)")

# %%% 5. Submit GBIF occurrence download ---------------------------------------
# occ_download() runs server-side on GBIF, so a single request handles all 
# species at once. If a previous run already completed the download, skip
# straight to loading the saved raw data (Step 6).
raw_rds <- file.path(out_dir, "GBIF_occurrences_raw.rds")

if (file.exists(raw_rds) && !FORCE_RERUN) {
  
  log_msg("Found existing raw download at ", raw_rds, " -- loading from disk (skipping GBIF download)")
  dat_raw <- readRDS(raw_rds)
  
} else {
  
  download_key <- occ_download(
    pred_in("taxonKey", gbif_keys),
    pred_in("country", countries),
    pred_gte("year", year_range[1]),
    pred_lte("year", year_range[2]),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_not(pred_in("basisOfRecord", c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))),  # exclude fossils and ex situ (zoo/garden) specimens
    format = "SIMPLE_CSV"
  )
  
  log_msg("Download submitted: ", download_key,
          " -- waiting for GBIF to prepare it (status printed every 60s)...")
  
  # status_ping controls how often occ_download_wait() polls and prints status,
  # so background progress is visible in run_log.txt while you're away.
  occ_download_wait(download_key, status_ping = 60, quiet = FALSE)
  
  dl_path <- occ_download_get(download_key, path = out_dir, overwrite = TRUE)
  dat_raw <- occ_download_import(dl_path)
  
  log_msg(nrow(dat_raw), " raw occurrence records downloaded")
  
  saveRDS(dat_raw, raw_rds)
  
  # Save the citation - cite this in any resulting publication
  writeLines(format(gbif_citation(dl_path)), file.path(out_dir, "gbif_citation.txt"))
}

# Build lookup: input name -> GBIF accepted name, by comparing the input
# species list against the species names that actually appear in the downloaded
# data. GBIF's download labels records with its accepted species name, so any
# input name that doesn't appear verbatim in dat_raw$species was a synonym.
gbif_species <- unique(dat_raw$species)
extra_species <- setdiff(gbif_species, species_list)   # names in download but not in input
missing_input <- setdiff(species_list, gbif_species)   # input names not in download (synonym or no records)

# Match each extra (accepted) name back to the input (synonym) name via the
# backbone's usageKey: both the input synonym and the accepted name share the
# same accepted usageKey, so we can join on that.
name_lookup <- tibble(input_name = character(), gbif_name = character())
if (length(extra_species) > 0 && length(missing_input) > 0) {
  # For each extra species, look up its usageKey, then find which input name
  # mapped to the same key
  extra_keys <- map_dfr(extra_species, function(sp) {
    bb <- name_backbone(name = sp, kingdom = "Plantae")
    tibble(gbif_name = sp, key = bb$usageKey)
  })
  input_keys <- accepted %>%
    transmute(input_name = verbatim_name, key = usageKey)
  name_lookup <- inner_join(input_keys, extra_keys, by = "key") %>%
    filter(input_name != gbif_name) %>%
    select(input_name, gbif_name) %>%
    distinct()
}

if (nrow(name_lookup) > 0) {
  log_msg(nrow(name_lookup), " input names are synonyms in GBIF's backbone:")
  for (i in seq_len(nrow(name_lookup))) {
    log_msg("  ", name_lookup$input_name[i], " -> ", name_lookup$gbif_name[i])
  }
}

# %%% 6. Prepare data for cleaning ----------------------------------------------
dat <- dat_raw %>%
  select(species, decimalLongitude, decimalLatitude, countryCode,
         individualCount, gbifID, family, taxonRank,
         coordinateUncertaintyInMeters, year, basisOfRecord,
         institutionCode, datasetKey) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

# ISO2 -> ISO3, required by cc_coun()
dat$countryCode <- countrycode(dat$countryCode, origin = "iso2c", destination = "iso3c")

dat <- as.data.frame(dat)
dat$row_id <- seq_len(nrow(dat))   # stable index, used to realign parallel chunk results

log_msg(nrow(dat), " records ready for cleaning")

# %%% 7. CoordinateCleaner record-level tests -----------------------------------
# cc_val/cc_equ/cc_cap/cc_cen/cc_coun/cc_sea/cc_zero/cc_dupl are single
# vectorized passes over the whole table (not per-species loops), so they run
# sequentially here.
log_msg("Running record-level CoordinateCleaner tests...")

clean <- dat %>%
  cc_val() %>%
  cc_equ() %>%
  cc_cap() %>%
  cc_cen() %>%
  cc_coun(iso3 = "countryCode") %>%
  cc_sea() %>%
  cc_zero() %>%
  cc_dupl()

log_msg(nrow(dat) - nrow(clean), " records removed by record-level tests; ", nrow(clean), " remaining")

# %%% 7b. cc_outl(), parallelized by species -------------------------------------
# cc_outl() evaluates each species' records independently (distance to that
# species' other points).
outl_rds <- file.path(out_dir, "cc_outl_flags.rds")

if (file.exists(outl_rds) && !FORCE_RERUN) {
  
  log_msg("Found existing cc_outl() flags at ", outl_rds, " -- loading from disk")
  outl_flags <- readRDS(outl_rds)
  
} else {
  
  log_msg("Running cc_outl() across ", length(unique(clean$species)),
          " species (", N_WORKERS, " workers)...")
  
  species_chunks <- split(unique(clean$species),
                          cut(seq_along(unique(clean$species)), N_WORKERS, labels = FALSE))
  
  outl_results <- future_map(species_chunks, function(sp_chunk) {
    sub <- clean[clean$species %in% sp_chunk, ]
    flags <- cc_outl(sub, species = "species", value = "flagged", verbose = FALSE)
    tibble(row_id = sub$row_id, outl_flag = flags)
  }, .progress = FALSE)
  
  outl_flags <- bind_rows(outl_results) %>% arrange(row_id)
  saveRDS(outl_flags, outl_rds)
}

clean <- clean %>%
  left_join(outl_flags, by = "row_id") %>%
  filter(outl_flag) %>%
  select(-outl_flag)

log_msg(nrow(clean), " records remaining after cc_outl()")

# %%% 8. Coordinate precision filter --------------------------------------------
# Records with no reported uncertainty are kept (NA), per the CoordinateCleaner
# tutorial's convention; recommended to inspect these separately.
clean <- clean %>%
  filter(is.na(coordinateUncertaintyInMeters) |
           (coordinateUncertaintyInMeters / 1000) <= coord_uncertainty_km)

log_msg(nrow(clean), " records remaining after precision filter (<= ",
        coord_uncertainty_km, " km)")

# %%% 9. Dataset-level tests, parallelized by dataset -----------------------------
# cd_ddmm() and cd_round() group by datasetKey, so chunks are built by
# splitting whole datasets across workers.
dataset_rds <- file.path(out_dir, "dataset_level_flags.rds")

if (file.exists(dataset_rds) && !FORCE_RERUN) {
  
  log_msg("Found existing dataset-level flags at ", dataset_rds, " -- loading from disk")
  saved <- readRDS(dataset_rds)
  ddmm_report  <- saved$ddmm_report
  round_report <- saved$round_report
  ddmm_flags   <- saved$ddmm_flags
  round_flags  <- saved$round_flags
  
} else {
  
  log_msg("Running cd_ddmm() and cd_round() across ", length(unique(clean$datasetKey)),
          " datasets (", N_WORKERS, " workers)...")
  
  dataset_ids <- unique(clean$datasetKey)
  dataset_chunks <- split(dataset_ids,
                          cut(seq_along(dataset_ids), N_WORKERS, labels = FALSE))
  
  chunk_results <- future_map(dataset_chunks, function(ds_chunk) {
    sub <- clean[clean$datasetKey %in% ds_chunk, ]
    
    # tryCatch: cd_ddmm() can fail on datasets with very few records (binom.test
    # receives invalid n < x). Treat such chunks as passing (all TRUE).
    ddmm_ds <- tryCatch({
      res <- cd_ddmm(sub, lon = "decimalLongitude", lat = "decimalLatitude",
                     ds = "datasetKey", diagnostic = TRUE, diff = 1, value = "dataset")
      tibble::rownames_to_column(res, var = "datasetKey")
    }, error = function(e) {
      message("  cd_ddmm() failed for a chunk (", length(ds_chunk), " datasets): ", conditionMessage(e))
      NULL
    })
    ddmm_fl <- tryCatch({
      cd_ddmm(sub, lon = "decimalLongitude", lat = "decimalLatitude",
              ds = "datasetKey", diff = 1, value = "flagged")
    }, error = function(e) rep(TRUE, nrow(sub)))
    
    round_ds <- tryCatch({
      cd_round(sub, lon = "decimalLongitude", lat = "decimalLatitude",
               ds = "datasetKey", T1 = 7, value = "dataset", graphs = FALSE)
    }, error = function(e) {
      message("  cd_round() failed for a chunk (", length(ds_chunk), " datasets): ", conditionMessage(e))
      NULL
    })
    round_fl <- tryCatch({
      cd_round(sub, lon = "decimalLongitude", lat = "decimalLatitude",
               ds = "datasetKey", T1 = 7, value = "flagged", graphs = FALSE)
    }, error = function(e) rep(TRUE, nrow(sub)))
    
    list(
      ddmm_report  = ddmm_ds,
      round_report = round_ds,
      flags        = tibble(row_id = sub$row_id, ddmm_flag = ddmm_fl, round_flag = round_fl)
    )
  }, .progress = FALSE)
  
  ddmm_report  <- bind_rows(map(chunk_results, "ddmm_report"))
  round_report <- bind_rows(map(chunk_results, "round_report"))
  flags_all    <- bind_rows(map(chunk_results, "flags")) %>% arrange(row_id)
  ddmm_flags   <- flags_all$ddmm_flag
  round_flags  <- flags_all$round_flag
  
  write_csv(ddmm_report, file.path(out_dir, "ddmm_conversion_flags_by_dataset.csv"))
  write_csv(round_report, file.path(out_dir, "rasterized_sampling_flags_by_dataset.csv"))
  
  saveRDS(list(ddmm_report = ddmm_report, round_report = round_report,
               ddmm_flags = ddmm_flags, round_flags = round_flags),
          dataset_rds)
}

# Remove records belonging to datasets flagged by either test.
# NOTE: these are dataset-level tests, so review ddmm_conversion_flags_by_dataset.csv
# and rasterized_sampling_flags_by_dataset.csv before treating this as final -
# a flagged dataset may still contain a mix of good and bad records.
clean_final <- clean[ddmm_flags & round_flags, ]

log_msg(nrow(clean) - nrow(clean_final),
        " records removed for belonging to a flagged dataset (ddmm/rasterized tests)")
log_msg(nrow(clean_final), " final cleaned records")

# %%% 10. Save outputs -----------------------------------------------------------
# Rename column names, into a format more useful for GFS
clean_final <- clean_final %>%
  select(-row_id) %>%
  rename(decimallongitude = decimalLongitude, decimallatitude = decimalLatitude)

species_dir <- file.path(out_dir, "species_csvs")
dir.create(species_dir, showWarnings = FALSE)

for (sp in unique(clean_final$species)) {
  sp_data <- clean_final[clean_final$species == sp, ]
  
  # If GBIF's accepted name differs from the input name, prepend the input name
  input_match <- name_lookup$input_name[name_lookup$gbif_name == sp]
  if (length(input_match) == 1 && input_match != sp) {
    file_prefix <- paste0(gsub(" ", "_", input_match), "_", gsub(" ", "_", sp))
  } else {
    file_prefix <- gsub(" ", "_", sp)
  }
  
  sp_file <- file.path(species_dir, paste0(file_prefix, "_", nrow(sp_data), "n_",
                                           format(Sys.Date(), "%Y-%m-%d"), ".csv"))
  write_csv(sp_data, sp_file)
}

elapsed <- round(difftime(Sys.time(), script_start, units = "mins"), 1)
log_msg("=== Done in ", elapsed, " minutes. Cleaned data written to: ",
        species_dir, " (one CSV per species) ===")

# Reconcile: which input species have no records in the final output, and why?
# Exclude synonym species that DO have output under their accepted name
synonym_inputs <- name_lookup$input_name[name_lookup$gbif_name %in% unique(clean_final$species)]
missing_spp <- setdiff(species_list, c(unique(clean_final$species), synonym_inputs))
if (length(missing_spp) > 0) {
  log_msg(length(missing_spp), " species from the input list have no final output:")
  missing_info <- tibble(species = missing_spp) %>%
    mutate(
      in_gbif_download  = species %in% unique(dat_raw$species),
      after_cc_records   = species %in% unique(dat$species),
      after_cc_pipeline = species %in% unique(clean$species),
      after_precision   = species %in% unique(clean_final$species)
    )
  for (i in seq_len(nrow(missing_info))) {
    row <- missing_info[i, ]
    if (!row$in_gbif_download) {
      reason <- "no records returned by GBIF download (may be a synonym; check synonym_mapping.csv)"
    } else if (!row$after_cc_records) {
      reason <- "all records lacked coordinates"
    } else if (!row$after_cc_pipeline) {
      reason <- "all records removed by CoordinateCleaner tests (cc_val/cc_equ/.../cc_outl)"
    } else {
      reason <- "all records removed by precision or dataset-level filters"
    }
    log_msg("  - ", row$species, ": ", reason)
  }
  write_csv(missing_info, file.path(out_dir, "species_missing_from_output.csv"))
} else {
  log_msg("All ", length(species_list), " input species have records in the final output")
}
