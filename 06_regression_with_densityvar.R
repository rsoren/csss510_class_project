#
# 06_regression_with_densityvar.R
#
# Reed Sorensen
# June 2016
#


rm(list = ls())

library(foreign)
library(spatstat)
library(SpatialEpi)
library(dplyr)



result1 <- readRDS("tmp3.RDS")

parameter_names_v2 <- c(
  "constant", "age", "birthorder", "married", "pca", "rural",
  "sigma_hosp", "sigma_nonhosp"
)

matrix(
  c(round(result1$par, digits = 4), round(sqrt(diag(solve(-1 * result1$hessian))), digits = 4)),
  dimnames = list(parameter_names_v2, c("Estimate", "Std. Error")),
  ncol = 2
)


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

# this version makes no distinctions by facility type or urbanicity
df_list <- list(
  df_alldat = df
)


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


## get surface var for final data set

param_sigmas_post <- c(0.7606, 0.8907)

df_cluster_post <- df_cluster

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
    x_var = df_cluster_post$xvar,
    y_var = df_cluster_post$yvar)

  tmp[is.na(tmp) | tmp < 0] <- 0

  # # rescale to original
  tmp * (max(df_list[[i]]$deliveryready, na.rm = TRUE) / max(tmp))
}))


df_cluster_post$surface_var <- apply(surface_dat_post, MARGIN = 1, max)


## get density variable for final data set





df_post <- left_join(df_individual, df_cluster_post, by = c("clusterid" = "cluster_id")) %>%
  filter(complete.cases(.))



