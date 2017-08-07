#
# 07_stacked_density_version.R
#
# Reed Sorensen
# June 2017
#


library(foreign)
library(dplyr)
library(data.table)
library(ggplot2)
library(SpatialEpi)
library(Formula)
library(KernSmooth)
library(lattice)
library(lme4)
library(clusterSEs)



rm(list = ls())


#####
# DATA PREP

# -- read in individual and cluster data

df_individual <- read.csv("data/haitibirthclean3.csv") %>%
  mutate(wealth = wealth / 10000)

df_cluster <- read.dbf("data/DHS GPS/HTGE61FL.dbf") %>%
  dplyr::select(cluster_id = DHSCLUST, latitude = LATNUM, longitude = LONGNUM) %>%
  filter(latitude != 0 & longitude != 0)

df_cluster[, c("xvar", "yvar")] <- latlong2grid( # convert lat/long to km
  df_cluster[, c("longitude", "latitude")]
)

# # check which individual cluster IDs aren't in the 'df_cluster' dataset
# unique(df_individual$clusterid[!df_individual$clusterid %in% df_cluster$cluster_id])
# # 13  14 179 297 319 338 400
#
# # vice versa
# unique(df_cluster$cluster_id[!df_cluster$cluster_id %in% df_individual$clusterid])
# # 18  29 110 122 164 223 232 267 295 298 412


# -- read in facility data

dat1 <- read.csv("data/facility_data_anc.csv")
dat2 <- read.csv("data/haitideliveryindicators.csv")

df <- left_join(dat1, dat2, by = c("facility_id" = "facil")) %>%
  filter(!is.na(deliveryready))

# convert lat/long to kilometers
df[, c("xvar", "yvar")] <- latlong2grid(df[, c("longitude", "latitude")])


# repeat rows of the facility data, where
#   the number of reps is the value of 'deliveryready'
df2 <- data.table::rbindlist(lapply(unique(df$facility_id), function(x) {
  tmp <- subset(df, facility_id == x)
  tmp[rep(1, times = as.integer(tmp$deliveryready)), ]
}))


# -- define extent and resolution for spatial analysis
pixel_length_km <- 0.5

xlo <- min(df$xvar) - pixel_length_km
xhi <- max(df$xvar) + pixel_length_km
ylo <- min(df$yvar) - pixel_length_km
yhi <- max(df$yvar) + pixel_length_km

n_xpoints <- ceiling((xhi - xlo) / pixel_length_km)
n_ypoints <- ceiling((yhi - ylo) / pixel_length_km)
n_points <- c(n_xpoints, n_ypoints)



#####
# FUNCTIONS

# function for getting cluster-specific estimates for a surface, given a bandwidth
# -- this is used for KDE with 'facilitybirth' values included (df2), and
#    KDE with density only (df)
# -- setting up the function this way allows for separate bandwidths by facility type/location

get_surface_estimate <- function(bw, dataset, by_cluster = TRUE) {

  # bw <- 0.8; use_deliveryready <- FALSE # dev variable
  require(akima); require(KernSmooth)

  # fit the surface
  surface1 <- KernSmooth::bkde2D(
    x = as.data.frame(dataset)[, c("xvar", "yvar")],
    bandwidth = bw,
    gridsize = n_points
  )

  if (by_cluster) {
    out <- akima::bilinear(
      x = surface1$x1, y = surface1$x2, z = surface1$fhat,
      x0 = df_cluster$xvar, y0 = df_cluster$yvar
    )

  } else { out <- surface1$fhat }

  return(out)

}




# function for getting a 'glm' object with logistic regression results
# -- optionally also returns the individual-level data used to fit the regression

get_model_fit <- function(param_vals, formula_var, return_individual_dat = FALSE) {

  # param_vals <- c(10,15) # dev
  # formula_var <- formula1 # dev
  # return_individual_dat = FALSE # dev
  # density_surface_option <- "subtract" # dev

  surface_dat <- lapply(1:length(facility_type_list), function(i) {
    x <- facility_type_list[[i]]
    df_orig <- subset(df, eval(x))
    df_repeated <- subset(df2, eval(x))
    density_and_value <- get_surface_estimate(param_vals[i], df_repeated)

    # three options:
    # do weighted_density minus nonweighted_density
    if (density_surface_option == "subtract") {
      density_only <- get_surface_estimate(param_vals[i], df_orig)
      out <- list(covar1 = density_and_value$z - density_only$z)

    # return only weighted_density; uses nonweighted_density as covariate later
    } else if (density_surface_option == "as_covariate") {
      out <- list(
        covar1 = density_and_value$z,
        covar2 = get_surface_estimate(param_vals[i], df_orig)$z
      )

    # return only weighted_density
    } else if (density_surface_option == "none") {
      out <- list(covar1 = density_and_value$z)
    }

    return(out)
  })

  surface_dat_covar1 <- do.call("cbind", lapply(surface_dat, function(x) x[["covar1"]]))
  df_cluster$surface_var <- apply(surface_dat_covar1, MARGIN = 1, max)

  # create another variable 'surface_dat_covar2' if "as_covariate" option is selected
  if (density_surface_option == "as_covariate") {
    surface_dat_covar2 <- do.call("cbind", lapply(surface_dat, function(x) x[["covar2"]]))
    df_cluster$density_var <- apply(surface_dat_covar2, MARGIN = 1, max)
  }

  df_individual2 <- df_cluster %>%
    right_join(df_individual, by = c("cluster_id" = "clusterid")) %>%
    mutate(svyweight = as.integer(svyweight / 10000) )

  out <- list( fit = glm(
    formula = formula_var,
    data = df_individual2,
    family = "binomial",
    weights = svyweight
  ))

  if (return_individual_dat) out <- append(out, list(individual_dat = df_individual2))

  return(out)

}

# function for viewing the estimate and standard error of parameters
#   after optim() finishes
print_result <- function(res, stval) {
  matrix(
    c(round(res$par, digits = 4), round(sqrt(diag(solve(-1 * res$hessian))), digits = 4)),
    dimnames = list(names(stval), c("Estimate", "Std. Error")),
    ncol = 2
  )
}



#####
# MODEL SPECIFICATIONS

# 1. Include 'rural' in individual-level logistic regression formula
#    -- don't estimate a separate bandwidth by urbanicity
#
# specification_name <- "2bw_age_married_bo_pca_rural_surface_subtractdensity"
# formula1 <- Formula(
#   facilitybirth ~ age + married + birthorder + pca + rural + surface_var)
# density_surface_option <- "subtract" # "none", "subtract", "as_covariate"
#
# facility_type_list <- list(
#   quote(hospital == 1),
#   quote(hospital == 0)
# )
#
# parameter_stval <- c(bw1 = 6, bw2 = 11)


# 2. Don't include 'rural' in individual-level logistic regression formula
#    -- estimate separate bandwidth by urbanicity (among clinics, not hospitals)
#
# specification_name <- "3bw_age_married_bo_pca_surface_subtractdensity"
# formula1 <- Formula(facilitybirth ~ age + married + birthorder + pca + surface_var)
# density_surface_option <- "subtract" # "none", "subtract", "as_covariate"
#
# facility_type_list <- list(
#   quote(hospital == 1), # not enough rural hospitals (18) to split hospitals by urbanicity
#   quote(hospital == 0 & rural == 1),
#   quote(hospital == 0 & rural == 0)
# )
# parameter_stval <- c(bw_h=8, bw_cr=35, bw_cu=11)


# 3. Include 'rural' in individual-level logistic regression formula
#    -- estimate separate bandwidth by urbanicity (among clinics, not hospitals)
#
# specification_name <- "3bw_age_married_bo_pca_rural_surface_subtractdensity"
# formula1 <- Formula(facilitybirth ~ age + married + birthorder + pca + rural + surface_var)
# density_surface_option <- "subtract" # "none", "subtract", "as_covariate"
#
# facility_type_list <- list(
#   quote(hospital == 1), # not enough rural hospitals (18) to split hospitals by urbanicity
#   quote(hospital == 0 & rural == 1),
#   quote(hospital == 0 & rural == 0)
# )
# parameter_stval <- c(bw_h=8, bw_cr=35, bw_cu=11) # values based on results of previous run


# 4. Estimate separate bandwidths for urban/rural only, not hospital/clinic
#    -- use all other covariates
#
specification_name <- "2bw_byurbanicity_age_married_bo_pca_rural_surface_subtractdensity_useweights"
formula1 <- Formula(facilitybirth ~
    age + married + birthorder + pca + rural + surface_var)
density_surface_option <- "subtract" # "none", "subtract", "as_covariate"

facility_type_list <- list(
  quote(rural == 1),
  quote(rural == 0)
)
parameter_stval <- c(bw_rural=8, bw_urban=8)



# 5. Make no distinctions by health facility or urbanicity
#    -- Include one covariate in the final regression:
#         readiness-weighted density minus plain density
#
# specification_name <- "1bw_age_married_bo_pca_rural_surface_subtractdensity"
# formula1 <- Formula(facilitybirth ~ age + married + birthorder + pca + rural + surface_var)
# density_surface_option <- "subtract" # "none", "subtract", "as_covariate"
#
# facility_type_list <- list(
#   quote(!is.na(rural))
# )
# parameter_stval <- c(bw_all=10)


# 6. Include 'rural' in individual-level logistic regression formula
#    -- estimate separate bandwidth by urbanicity (among clinics, not hospitals)
#    -- include KDE density (non-weighted) as a covariate, not subtraction
# specification_name <- "2bw_byurbanicity_age_married_bo_pca_rural_surface_densitycovar_useweights"
# formula1 <- Formula(facilitybirth ~
#     age + married + birthorder + pca + rural + surface_var + density_var)
# density_surface_option <- "as_covariate" # "none", "subtract", "as_covariate"
#
# facility_type_list <- list(
#   quote(rural == 1),
#   quote(rural == 0)
# )
# parameter_stval <- c(bw_rural=8, bw_urban=8)


# 7. pca and spatial vars are only covariates in logistic regression
#    -- estimate separate bandwidth by urbanicity (among clinics, not hospitals)
#    -- include KDE density (non-weighted) as a covariate, not subtraction
#
# specification_name <- "2bw_byurbanicity_pca_surface_densitycovar"
# formula1 <- Formula(facilitybirth ~ pca + surface_var + density_var)
# density_surface_option <- "as_covariate" # "none", "subtract", "as_covariate"
#
# facility_type_list <- list(
#   quote(rural == 1),
#   quote(rural == 0)
# )
# parameter_stval <- c(bw_rural=8, bw_urban=35)





#####
# RUN MODEL

llk_logit <- function(param) {
  # param <- c(5,5,5,5) # dev variable
  if (any(param <= 0.01)) return(-100000)
  fit1 <- get_model_fit(param_vals = param, formula_var = formula1)
  as.numeric(strsplit(as.character(logLik(fit1[["fit"]])), split = ' ')[[1]]) # extract likelihood
}

system.time(result1 <- optim(
  par = parameter_stval, fn = llk_logit,
  method = "BFGS", hessian = T, control = list(fnscale = -1)
))


#####
# RESULTS

load_previous_result <- TRUE

if (load_previous_result) {
  result1 <- readRDS(
    "results/2bw_byurbanicity_age_married_bo_pca_rural_surface_subtractdensity_useweights.RDS" )[[2]]
}

(displayed_result <- print_result(res = result1, stval = parameter_stval))


# re-run using the estimated sigmas to get final logistic regression results

fit2 <- get_model_fit(
  param_vals = result1$par,
  formula_var = formula1,
  return_individual_dat = TRUE
)

summary(fit2[["fit"]])


# get results using standard errors that account for clustering
fit3 <- clusterSEs::cluster.bs.glm(
  mod = fit2[["fit"]],
  dat = fit2[["individual_dat"]],
  cluster = ~ cluster_id,
  report = TRUE
)


# write result to disk
saveRDS(
  object = list(fit2, result1, fit3, displayed_result),
  file = paste0("results/", specification_name, ".RDS")
)



#####
# MAP

tmp <- readRDS("results/3bw_age_married_bo_pca_rural_surface_subtractdensity.RDS")
tmp_glm_fit <- tmp[[1]][["fit"]]
tmp_dat <- tmp[[1]][["individual_dat"]]
tmp_optim_object <- tmp[[2]]
tmp_bw_results <- tmp[[3]]


# function for returning data for the entire spatial grid, in the form of a matrix
# -- not just cluster-specific values as with the function 'get_surface_estimates'

surface_matrices <- function(bw, subtract_density) {

  do.call("cbind", lapply(1:length(facility_type_list), function(i) {

    x <- facility_type_list[[i]]
    df_orig <- subset(df, eval(x))
    df_repeated <- subset(df2, eval(x))
    density_and_value <- as.vector(
      get_surface_estimate(final_bw[i], df_repeated, by_cluster = FALSE) )

    if (subtract_density) {
      density_only <- as.vector(
        get_surface_estimate(final_bw[i], df_orig, by_cluster = FALSE) )
      out <- density_and_value - density_only
    } else if (!subtract_density) {
      out <- density_and_value
    }

    return(out)
  }))
}

# take the max value at each pixel, across the 3 surfaces
combined_surface <- apply(
  X = surface_matrices(bw = tmp_optim_object$par, subtract_density = TRUE),
  MARGIN = 1, max
)

# reshape spatial data into a matrix
combined_surface_matrix <- t(matrix(combined_surface, ncol = n_xpoints, byrow = TRUE))
combined_surface_matrix[combined_surface_matrix < 0] <- 0

# get the outline of haiti, and change the coordinates to match the map
haiti <- readRDS("data/HTI_adm0.rds")


# map the country
corners <- as.matrix(expand.grid(
  x = c(xlo, xhi),
  y = c(ylo, yhi) )) %>%
  grid2latlong(.)

library(maps)
library(mapdata)



pal <- colorRampPalette(c("white", "navyblue"), space = "rgb")
levelplot(combined_surface_matrix, xlab="", ylab="",
  row.values = seq(min(corners$x), max(corners$x), length.out = n_xpoints),
  column.values = seq(min(corners$y), max(corners$y), length.out = n_ypoints),
  col.regions = pal(20), add = TRUE) +
  layer(sp.polygons(out))

maps::map("world", "haiti", plot = TRUE)



# Instead of using lattice::levelplot directly, you can use the raster
# package for your regular grid data + projection information, and the
# rasterVis package rasterVis::levelplot.   You can then either use
# panel.spplot to overlay your state outline polygons, or use
# latticeExtra:layer and spplot() to add the state outlines.
# You may need to use spTransform() on your state outline polygons or
# projectRaster on the grid to get them both in the same desired projection.

library(rasterVis)
library(sp)

# Download States boundaries (might take time)
out <- getData('GADM', country='Haiti', level=1)

# Extract California state
California <- out[out$NAME_1 %in% 'California',]

# Plot raster and California:
levelplot(RAD2012.all) +
  layer(sp.polygons(California))


# map a zoomed in area
levelplot(combined_surface_matrix, xlab="", ylab="",
  col.regions = pal(20),
  xlim = c(300, 470),
  ylim = c(200, 350)
)



