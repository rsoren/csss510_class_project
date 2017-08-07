#
# scratch.R
#
# Reed Sorensen
# November 2016
#

dhs_path <- "C:/Users/rsoren/Documents/prog/repos/csss510_class_project/data/dhs_2012_haiti/"
spa_path <- "C:/Users/rsoren/Documents/prog/repos/csss510_class_project/data/spa_2013_haiti/"

require(foreign)
# HTFC6BFLSP

dhs_files <- list.files(paste0(dhs_data), pattern = "DTA", recursive = TRUE)
spa_files <- list.files(paste0(spa_data), pattern = "DTA", recursive = TRUE)

dhs_dat <- lapply(dhs_files, function(x) read.dta(paste0(dhs_path, x)))
spa_dat <- lapply(spa_files, function(x) read.dta(paste0(spa_path, x)))

dhs_names <- do.call("c", lapply(dhs_dat, names))
spa_names <- do.call("c", lapply(spa_dat, names))

dhs_names[dhs_names == "LATNUM"]
spa_names[spa_names == "latnum"]

require(haven)
df <- haven::read_dta("data/HT13AFC_Haiti_SPA_FacilityInventory.dta")
df <- read.dbf()



