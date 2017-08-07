#
# 02_mle_example.R
#
# Reed Sorensen
# November 2016
#

rm(list = ls())

library(foreign)
library(spatstat)
library(SpatialEpi)
library(dplyr)


#####
# read in individual and cluster data

df_individual <- read.csv("data/haitibirthclean3.csv") %>%
  mutate(wealth = wealth / 10000)

df_cluster <- read.dbf("data/DHS GPS/HTGE61FL.dbf") %>%
  dplyr::select(cluster_id = DHSCLUST, latitude = LATNUM, longitude = LONGNUM) %>%
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

# df_list <- list(
#   df_hosp_urban = subset(df, hospital == 1 & rural == 0),
#   df_hosp_rural = subset(df, hospital == 1 & rural == 1),
#   df_nonhosp_urban = subset(df, hospital == 0 & rural == 0),
#   df_nonhosp_rural = subset(df, hospital == 0 & rural == 1)
# )

# # this version only differentiates between hospital and non-hospital
df_list <- list(
  df_hosp = subset(df, hospital == 1),
  df_nonhosp = subset(df, hospital == 0)
)

# this version makes no distinctions by facility type or urbanicity
# df_list <- list(
#   df_alldat = df
# )


#####
# set pixel and grid size

pixel_length <- 0.5

xlo <- min(df$xvar) - pixel_length
xhi <- max(df$xvar) + pixel_length
ylo <- min(df$yvar) - pixel_length
yhi <- max(df$yvar) + pixel_length

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
    dplyr::select(xvar, yvar, deliveryready)

  # tmp$deliveryready <- 67 # counterfactual; max observed in dataset

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

# # version with 4 sigmas for facility type and urbanicity; wealth/educ PCA
# parameter_names <- c(
#   "constant", "wealth_education_PCA", "surface_var",
#   "sigma_hosp_urban", "sigma_hosp_rural", "sigma_nonhosp_urban", "sigma_nonhosp_rural"
# )
# parameter_stval <- c(-1.0043, 1.0019, 0.1063, 1.5, 3.0, 0.35, 0.7)


# PCA, age, surface_var and 1 KDE sigma
parameter_names <- c(
  "constant", "age", "birthorder", "married", "pca", "rural",
  # "sigma_hosp_urban", "sigma_hosp_rural", "sigma_clinic_urban", "sigma_clinic_rural"
  "sigma_hosp", "sigma_nonhosp"
  # "sigma_hosp_rural", "sigma_nonhosp_urban", "sigma_nonhosp_rural"
  # "sigma1"
)

# parameter_stval <- c(-4.2, 0.08, -0.8, 0.18, 1.26, 2.13, 0.7, 0.7, 0.7, 0.7)
parameter_stval <- c(-4.2, 0.08, -0.8, 0.18, 1.26, 2.13, 0.7, 0.7)
# parameter_stval <- c(-4.2, 0.08, -0.8, 0.18, 1.26, 2.13, 0.7)

# # version with 2 sigmas for facility type; birthorder, age, married, wealth/educ pca
# parameter_names <- c(
#   "constant", "birthorder", "age", "married", "wealth_education_PCA", "surface_var",
#   "sigma_hosp", "sigma_nonhosp"
# )
#
# parameter_stval <- c(-1.2591, -0.2200, 0.0325, 0.0157,
#   0.8903, 0.0119, 1.0975, 0.3562)

length(parameter_names) == length(parameter_stval)


llk.logit <- function(param) {

  # param <- parameter_stval # dev variable

  n_sigmas <- 2 # set this depending on the model
  param_regression <- param[1:(length(parameter_stval)-n_sigmas)]
  param_sigmas <- param[(length(param_regression)+1):(length(parameter_stval))]

  if (any(param_sigmas <= 0.01)) return(-100000)

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

    # sigma_var_tmp <- ifelse(i==1, 0.5, param_sigmas[(i-1)]) # changed from 1 to 0.5

    tmp <- get_surface_points(
      surface = df_list2[[i]],
      sigma_var = param_sigmas[i],
      # sigma_var = sigma_var_tmp,
      x_var = df_cluster_tmp$xvar,
      y_var = df_cluster_tmp$yvar)

    tmp[is.na(tmp) | tmp < 0] <- 0

    # # rescale to original
    tmp * (max(df_list[[i]]$deliveryready, na.rm = TRUE) / max(tmp))
  }))


  df_cluster_tmp$surface_var <- apply(surface_dat, MARGIN = 1, sum) # this doesn't rescale

  dat2 <- left_join(df_individual, df_cluster_tmp, by = c("clusterid" = "cluster_id")) %>%
    filter(complete.cases(.)) # figure out why 25 out of 1991 are NA; they aren't anymore

  y <- dat2$facilitybirth
  # x <- cbind(dat2$pca, dat2$surface_var)
  x <- cbind(dat2$age, dat2$birthorder, dat2$married,
    dat2$pca, dat2$rural, dat2$surface_var)

  os <- rep(1,length(x[,1]))
  x2 <- cbind(os,x)
  # b <- param_regression[1:ncol(x2)]
  b <- c(param_regression, 1)
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

# saveRDS(result1, "tmp.RDS") # version with baseline sigma = 1 (urban hospital)
# saveRDS(result1, "tmp2.RDS") # version with baseline sigma = 0.5 (urban hospital)
# saveRDS(result1, "tmp3.RDS") # version with 2 sigma, before adding density covariate
# saveRDS(result1, "tmp4.RDS") # version with 4 sigma, before adding density covariate
# saveRDS(result1, "tmp5.RDS") # version with 1 sigma, before adding density covariate
saveRDS(result1, "tmp6.RDS") # version with 2 sigma, before adding density, summing surfaces


result1 <- readRDS("tmp3.RDS")




# tmp <- readRDS("data/results1.RDS")
# tmp2 <- readRDS("data/results_hosp_nonhosp.RDS") # covars: wealth and educat
# tmp3 <- readRDS("data/results_onesigma.RDS") # covars: wealth and educat
# tmp4 <- readRDS("data/results_hosp_nonhosp_wealthonly.RDS") # covars: wealth
# tmp5 <- readRDS("data/results_hosp_nonhosp_educonly.RDS") # covars: educ
# tmp6 <- readRDS("data/results_2sigmas_wealth_education.RDS")
# tmp7 <- readRDS("data/results_2sigmas_pca_v2.RDS")
# tmp8 <- readRDS("data/results_4sigmas_pca.RDS")
# tmp9 <- readRDS("data/results_2sigmas_pca_and_all_others_v2.RDS")
# tmp10 <- readRDS("data/new_results_2sigmas_age_birthorder_pca.RDS")
tmp11 <- readRDS("data/new_results_4sigmas_age_birthorder_pca.RDS")
#
result1 <- tmp11
#
saveRDS(result1, "data/new_results_4sigmas_age_birthorder_pca.RDS")




## regression


param_sigmas_post <- c(0.7606, 0.8907)


get_surface_points <- function(surface, sigma_var, x_var, y_var) {
  surface1 <- Smooth(surface, sigma = sigma_var)
  get_pixel_value1 <- spatstat::as.function.im(surface1)
  # get pixel values (and fix where the function returns an empty vector)
  tmp <- mapply(get_pixel_value1, x = x_var, y = y_var)
  tmp[unlist(lapply(tmp, function(x) length(x)==0))] <- NA
  return(unlist(tmp))
}

surface_dat_post <- do.call("cbind", lapply(1:length(df_list2), function(i) {

  # sigma_var_tmp <- ifelse(i==1, 0.5, param_sigmas[(i-1)]) # changed from 1 to 0.5

  tmp <- get_surface_points(
    surface = df_list2[[i]],
    sigma_var = param_sigmas_post[i],
    # sigma_var = sigma_var_tmp,
    x_var = df_cluster_tmp$xvar,
    y_var = df_cluster_tmp$yvar)

  tmp[is.na(tmp) | tmp < 0] <- 0

  # # rescale to original
  tmp * (max(df_list[[i]]$deliveryready, na.rm = TRUE) / max(tmp))
}))


df_cluster_post$surface_var <- apply(surface_dat_post, MARGIN = 1, max)


df_post <- left_join(df_individual, df_cluster_post, by = c("clusterid" = "cluster_id")) %>%
  filter(complete.cases(.))








