#
# 02_mle_example.R
#
# Reed Sorensen
# November 2016
#

rm(list = ls())

library(dplyr)
library(foreign)
library(spatstat)
library(SpatialEpi)


#####
# read in individual and cluster data

df_individual <- read.csv("data/haitibirthclean2.csv") %>%
  mutate(wealth = wealth / 10000)

df_cluster <- read.dbf("data/DHS GPS/HTGE61FL.dbf") %>%
  select(cluster_id = DHSCLUST, latitude = LATNUM, longitude = LONGNUM) %>%
  filter(latitude != 0 & longitude != 0)

df_cluster[, c("xvar", "yvar")] <- latlong2grid( # convert lat/long to km
  df_cluster[, c("longitude", "latitude")]
)


#####
# read in facility data

dat1 <- read.csv("data/facility_data_anc.csv")
dat2 <- read.csv("data/haitideliveryindicators.csv")

df <- left_join(dat1, dat2, by = c("facility_id" = "facil")) %>%
  filter(!is.na(deliveryready))

df[, c("xvar", "yvar")] <- latlong2grid(df[, c("longitude", "latitude")])


# split data by facility type and urbanicity
#   in order to estimate KDE sigma for each one

# # NOTE: there are only 18 rural hospitals, and
# # the model failed to give an estimate for sigma_hosp_rural

df_list <- list(
  df_hosp_urban = subset(df, hospital == 1 & rural == 0),
  df_hosp_rural = subset(df, hospital == 1 & rural == 1),
  df_nonhosp_urban = subset(df, hospital == 0 & rural == 0),
  df_nonhosp_rural = subset(df, hospital == 0 & rural == 1)
)

# # this version only differentiates between hospital and non-hospital
# df_list <- list(
#   df_hosp = subset(df, hospital == 1),
#   df_nonhosp = subset(df, hospital == 0)
# )

# this version makes no distinctions by facility type or urbanicity
# df_list <- list(
#   df_alldat = df
# )


#####
# set pixel and grid size

pixel_length <- 0.5

zoom <- FALSE

if (zoom) {

  xlo <- min(df$xvar) + 100 - pixel_length
  xhi <- min(df$xvar) + 120 + pixel_length
  ylo <- min(df$yvar) + 30 - pixel_length
  yhi <- min(df$yvar) + 50 + pixel_length

} else {

  xlo <- min(df$xvar) - pixel_length
  xhi <- max(df$xvar) + pixel_length
  ylo <- min(df$yvar) - pixel_length
  yhi <- max(df$yvar) + pixel_length

}

x_points <- seq(xlo, xhi, by = pixel_length)
y_points <- seq(ylo, yhi, by = pixel_length)

window1 <- owin(
  xrange = c(min(x_points), max(x_points)),
  yrange = c(min(y_points), max(y_points)))

dat_grid <- expand.grid( # grid of zeros
  xvar = x_points,
  yvar = y_points,
  deliveryready = 0
)


#####
# process facility-level data

df_list2 <- lapply(df_list, function(x) {

  # x <- df_list[[1]] # dev variable

  tmp <- x %>% # keep points within (potentially zoomed) area
    filter(xvar >= xlo & xvar <= xhi) %>%
    filter(yvar >= ylo & yvar <= yhi) %>%
    select(xvar, yvar, deliveryready)

  tmp2 <- rbind(tmp, dat_grid) # add grid of zeros

  tmp3 <- ppp( # create spatial data frame
    x = tmp2$xvar,
    y = tmp2$yvar,
    window = window1
  )

  marks(tmp3) <- tmp2$deliveryready # add facility values

  return(tmp3)

})


#####
# maximum likelihood function

# # version with hospital/non-hospital and urban/rural
# parameter_names <- c(
#   "constant", "wealth", "education", "surface_var",
#   "sigma_hosp_urban", "sigma_hosp_rural", "sigma_nonhosp_urban", "sigma_nonhosp_rural"
# )
# parameter_stval <- c(-1.3181, 0.0933, 0.7176, 0.0687, 1.6, 3.2, 0.377, 1.0)

# version with hospital/non-hospital distinction only
# parameter_names <- c(
#   "constant", "wealth", "education", "surface_var",
#   "sigma_hosp", "sigma_nonhosp"
# )
#
# parameter_stval <- c(-1.3181, 0.0993, 0.7176, 0.0687, 1.6034, 0.3774)


# # # version with hospital/non-hospital distinction only (wealth covariate only)
# parameter_names <- c(
#   "constant", "wealth", "surface_var",
#   "sigma_hosp", "sigma_nonhosp"
# )
#
# parameter_stval <- c(-0.39, 0.126, 0.021, 20, 5)

# # version with hospital/non-hospital distinction only (educat covariate only)
# parameter_names <- c(
#   "constant", "education", "surface_var",
#   "sigma_hosp", "sigma_nonhosp"
# )
#
# parameter_stval <- c(-2.1553, 1.1241, 0.3443, 1.8154, 0.7823)


# # version with no sigma distinctions
# parameter_names <- c(
#   "constant", "wealth", "education", "surface_var",
#   "sigma"
# )
#
# parameter_stval <- c(-0.39, 0.126, 0.7, 0.021, 5)





llk.logit <- function(param) {

  # param <- parameter_stval # dev variable

  # NOTE: change these depending on what's in the model
  param_regression <- param[1:4]
  param_sigmas <- param[5:8]

  if (any(param_sigmas <= 0.001)) return(-100000)

  df_cluster_tmp <- df_cluster

  get_surface_points <- function(surface, sigma_var, x_var, y_var) {
    surface1 <- Smooth(surface, sigma = sigma_var)
    get_pixel_value1 <- spatstat::as.function.im(surface1)
    # get pixel values (and fix where the function returns an empty vector)
    tmp <- mapply(get_pixel_value1, x = x_var, y = y_var)
    tmp[unlist(lapply(tmp, function(x) length(x)==0))] <- NA
    return(unlist(tmp))
  }

  surface_dat <- do.call("cbind", lapply(1:length(df_list2), function(i) {

    get_surface_points(
      surface = df_list2[[i]],
      sigma_var = param_sigmas[i],
      x_var = df_cluster_tmp$xvar,
      y_var = df_cluster_tmp$yvar )
  }))

  surface_dat[is.na(surface_dat) | surface_dat < 0] <- 0
  df_cluster_tmp$surface_var <- apply(surface_dat, MARGIN = 1, max)

  dat2 <- left_join(df_individual, df_cluster_tmp, by = c("clusterid" = "cluster_id")) %>%
    filter(complete.cases(.)) # figure out why 25 out of 1991 are NA

  y <- dat2$facilitybirth
  x <- cbind(dat2$wealth, dat2$educat, dat2$surface_var)
  # x <- cbind(dat2$educat, dat2$surface_var)

  os <- rep(1,length(x[,1]))
  x2 <- cbind(os,x)
  b <- param_regression[1:ncol(x2)]
  xb <- x2 %*% b
  sum(-1 * (y*log(1+exp(-xb)) + (1-y)*log(1+exp(xb))))
}



system.time(result1 <- optim(
  par = parameter_stval, fn = llk.logit,
  method = "BFGS", hessian = T, control = list(fnscale = -1)
))

matrix(
  c(round(result1$par, digits = 4), round(sqrt(diag(solve(-1 * result1$hessian))), digits = 4)),
  dimnames = list(parameter_names, c("Estimate", "Std. Error")),
  ncol = 2
)



# tmp <- readRDS("data/results1.RDS")
# tmp2 <- readRDS("data/results_hosp_nonhosp.RDS") # covars: wealth and educat
# tmp3 <- readRDS("data/results_onesigma.RDS") # covars: wealth and educat
# tmp4 <- readRDS("data/results_hosp_nonhosp_wealthonly.RDS") # covars: wealth
# tmp5 <- readRDS("data/results_hosp_nonhosp_educonly.RDS") # covars: educ
# tmp6 <- readRDS("data/results_facilitytype_wealtheducation.RDS")

result1 <- tmp6



# saveRDS(result1, "data/results_facilitytype_wealtheducation.RDS")






