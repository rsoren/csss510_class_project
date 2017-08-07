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
library(Formula)

library(KernSmooth)
library(plotrix)
library(heatmap3)
library(lattice)



rm(list = ls())

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

# convert lat/long to kilometers
df[, c("xvar", "yvar")] <- latlong2grid(df[, c("longitude", "latitude")])


# repeat rows of the facility data, where
#   the number of reps is the value of 'deliveryready'
df2 <- data.table::rbindlist(lapply(unique(df$facility_id), function(x) {
  tmp <- subset(df, facility_id == x)
  tmp[rep(1, times = as.integer(tmp$deliveryready)), ]
}))


pixel_length_km <- 0.5

xlo <- min(df$xvar) - pixel_length_km
xhi <- max(df$xvar) + pixel_length_km
ylo <- min(df$yvar) - pixel_length_km
yhi <- max(df$yvar) + pixel_length_km

n_xpoints <- ceiling((xhi - xlo) / pixel_length_km)
n_ypoints <- ceiling((yhi - ylo) / pixel_length_km)
n_points <- c(n_xpoints, n_ypoints)


get_surface_estimate <- function(bw, use_deliveryready = FALSE) {

  # bw <- 0.8; use_deliveryready <- FALSE # dev variable
  require(akima)

  tmp_dat <- df
  if (use_deliveryready) tmp_dat <- df2

  # fit the surface
  surface1 <- KernSmooth::bkde2D(
    x = as.data.frame(tmp_dat)[, c("xvar", "yvar")],
    bandwidth = bw,
    gridsize = n_points
  )

  akima::bilinear(
    x = surface1$x1,
    y = surface1$x2,
    z = surface1$fhat,
    x0 = df_cluster$xvar,
    y0 = df_cluster$yvar
  )
}




formula1 <- Formula(facilitybirth ~ pca + age + birthorder + rural + surface_var)

llk_logit <- function(param) {

  # param <- 0.8 # dev variable
  bandwidth_var <- param[1]

  density_and_value <- get_surface_estimate(bw = bandwidth_var, use_deliveryready = TRUE)
  density_only <- get_surface_estimate(bw = bandwidth_var)

  surface_var <- density_and_value$z - density_only$z

  df_cluster2 <- df_cluster %>%
    mutate(
      surface_var = surface_var ) %>%
    #   surface_var = rescale(surface_est$z, c(0, max(df$deliveryready))) ) %>%
    dplyr::select(cluster_id, surface_var)

  df_individual2 <- df_individual %>%
    left_join(df_cluster2, by = c("clusterid" = "cluster_id"))

  fit1 <- glm(
    formula = formula1,
    data = df_individual2,
    family = binomial(link = "logit")
  )

  # get log likelihood
  as.numeric(strsplit(as.character(logLik(fit1)), split = ' ')[[1]])

}


parameter_stval <- c(sigma = 0.8)

system.time(result1 <- optim(
  par = parameter_stval, fn = llk_logit,
  method = "BFGS", hessian = T, control = list(fnscale = -1)
))


matrix(
  c(round(result1$par, digits = 4), round(sqrt(diag(solve(-1 * result1$hessian))), digits = 4)),
  dimnames = list(names(parameter_stval), c("Estimate", "Std. Error")),
  ncol = 2
)


# re-run with the estimated sigma to get final coefficient estimates

density_and_value <- get_surface_estimate(bw = result1$par, use_deliveryready = TRUE)
density_only <- get_surface_estimate(bw = result1$par)

surface_var <- density_and_value$z - density_only$z

df_cluster2 <- df_cluster %>%
  mutate(
    surface_var = surface_var ) %>%
  #   surface_var = rescale(surface_est$z, c(0, max(df$deliveryready))) ) %>%
  dplyr::select(cluster_id, surface_var)

df_individual2 <- df_individual %>%
  left_join(df_cluster2, by = c("clusterid" = "cluster_id"))

fit2 <- glm(
  formula = formula1,
  data = df_individual2,
  family = binomial(link = "logit")
)

summary(fit2)


# # check which individual cluster IDs aren't in the 'df_cluster' dataset
# unique(df_individual$clusterid[!df_individual$clusterid %in% df_cluster$cluster_id])
# # 13  14 179 297 319 338 400
#
# # vice versa
# unique(df_cluster$cluster_id[!df_cluster$cluster_id %in% df_individual$clusterid])
# # 18  29 110 122 164 223 232 267 295 298 412



pal <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(surface1$fhat, xlab="", ylab="",
  col.regions = pal(10))

levelplot(surface1$fhat, xlab="", ylab="",
  col.regions = pal(10),
  xlim = c(170, 200),
  ylim = c(50, 70)
)







#-- get facility-readiness surface given a sigma value





#####

# # differentiate between hospital and non-hospital
# df_list <- list(
#   df_hosp = subset(df, hospital == 1),
#   df_nonhosp = subset(df, hospital == 0)
# )

# # no distinctions by facility type or urbanicity
df_list <- list(
  df_alldat = df
)













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













