# Rebuild the two real-data fixtures used by equivalence.R:
#   benchmarks/data/meck.gpkg       Mecklenburg County NC, ACS 2018-2022 tracts (293)
#   benchmarks/data/vancouver.gpkg  Vancouver, Burnaby and New Westminster BC, 2021 Census tracts (185)
# Both files are committed, so this script is only needed to regenerate them.
# Mecklenburg needs a Census API key (tidycensus::census_api_key()); Vancouver
# needs a CensusMapper key in CANCENSUS_API_KEY or CM_API_KEY (cancensus).
suppressPackageStartupMessages({library(sf); library(dplyr)})
out_dir <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "data")
dir.create(out_dir, showWarnings = FALSE)

# ---- Mecklenburg County, NC: median income, % bachelor's, median rent, median age
library(tidycensus)
options(tigris_use_cache = TRUE)
vars <- c(medinc = "B19013_001", pop = "B01003_001", bach = "B15003_022",
          edu_tot = "B15003_001", rent = "B25064_001", age = "B01002_001")
meck <- get_acs("tract", variables = vars, state = "NC", county = "Mecklenburg", year = 2022,
                geometry = TRUE, output = "wide", progress_bar = FALSE) |>
  transmute(GEOID, name = NAME, pop = popE, medinc = medincE,
            pct_bach = 100 * bachE / pmax(edu_totE, 1), rent = rentE, age = ageE) |>
  filter(pop > 0, !is.na(medinc), !is.na(rent), !is.na(age)) |>
  st_transform(32617) |>                      # UTM 17N
  filter(!st_is_empty(geometry))
st_write(meck, file.path(out_dir, "meck.gpkg"), delete_dsn = TRUE, quiet = TRUE)
cat("Mecklenburg County NC:", nrow(meck), "tracts\n")

# ---- Vancouver / Burnaby / New Westminster, BC: occupation shares (NOC 2021 broad categories)
library(cancensus)
key <- Sys.getenv("CANCENSUS_API_KEY"); if (!nzchar(key)) key <- Sys.getenv("CM_API_KEY")
options(cancensus.api_key = key, cancensus.cache_path = tempdir())
occ <- c(mgmt = "v_CA21_6570", business = "v_CA21_6573", sciences = "v_CA21_6576", health = "v_CA21_6579",
         education = "v_CA21_6582", arts = "v_CA21_6585", sales = "v_CA21_6588", trades = "v_CA21_6591",
         resources = "v_CA21_6594", manufacturing = "v_CA21_6597")
van <- get_census("CA21", regions = list(CSD = c("5915022", "5915025", "5915029")), level = "CT",
                  vectors = c(labour = "v_CA21_6567", occ), geo_format = "sf", quiet = TRUE) |>
  rename(pop = Population) |>
  filter(pop > 0, !is.na(labour), labour > 0)
for (a in names(occ)) van[[a]] <- round(100 * van[[a]] / van$labour, 2)
van <- van |> mutate(GEOID = GeoUID, name = paste("CT", GeoUID)) |>
  select(GEOID, name, pop, labour, all_of(names(occ))) |>
  st_transform(26910) |>                      # UTM 10N (NAD83)
  filter(!st_is_empty(geometry))
st_write(van, file.path(out_dir, "vancouver.gpkg"), delete_dsn = TRUE, quiet = TRUE)
cat("Vancouver/Burnaby/New Westminster:", nrow(van), "tracts (one island tract; equivalence.R keeps the main component)\n")
