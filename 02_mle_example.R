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

# facility data
dat1 <- read.csv("data/facility_data_anc.csv")
dat2 <- read.csv("data/haitideliveryindicators.csv")

df <- left_join(dat1, dat2, by = c("facility_id" = "facil")) %>%
  filter(!is.na(deliveryready))

df[, c("xvar", "yvar")] <- latlong2grid(df[, c("longitude", "latitude")])

df_hospital <- subset(df, hospital == 1)
df_nonhospital <- subset(df, hospital == 0)

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

df2_hospital <- df_hospital %>%
  filter(xvar >= xlo & xvar <= xhi) %>%
  filter(yvar >= ylo & yvar <= yhi) %>%
  select(xvar, yvar, deliveryready)

df2_nonhospital <- df_nonhospital %>%
  filter(xvar >= xlo & xvar <= xhi) %>%
  filter(yvar >= ylo & yvar <= yhi) %>%
  select(xvar, yvar, deliveryready)

# add grid of zeros
dat_grid <- expand.grid(xvar = x_points, yvar = y_points, deliveryready = 0)

df3_hospital <- rbind(df2_hospital, dat_grid)
df3_nonhospital <- rbind(df2_nonhospital, dat_grid)

window1 <- owin(
  xrange = c(min(x_points), max(x_points)),
  yrange = c(min(y_points), max(y_points)))

df_spatial_hospital <- ppp(
  x = df3_hospital$xvar,
  y = df3_hospital$yvar,
  window = window1
)

df_spatial_nonhospital <- ppp(
  x = df3_nonhospital$xvar,
  y = df3_nonhospital$yvar,
  window = window1
)

marks(df_spatial_hospital) <- df3_hospital$deliveryready
marks(df_spatial_nonhospital) <- df3_nonhospital$deliveryready

# individual and cluster data
df_dhs <- read.csv("data/haitibirthclean2.csv") %>%
  mutate(wealth = wealth / 10000)

df_cluster <- read.dbf("data/DHS GPS/HTGE61FL.dbf") %>%
  select(cluster_id = DHSCLUST, latitude = LATNUM, longitude = LONGNUM) %>%
  filter(latitude != 0 & longitude != 0)

df_cluster[, c("xvar", "yvar")] <- latlong2grid(
  df_cluster[, c("longitude", "latitude")]
)


# maximum likelihood function

llk.logit <- function(param) {

  param <- c(-0.39, 0.126, 2, 0.021, 40.0, 0.355) # dev variable
  # constant, wealth, education, surface_var, sigma1, sigma2

  param_regression <- param[1:4]
  param_sigma1 <- param[5]
  param_sigma2 <- param[6]

  if (param_sigma1 <= 0.001 | param_sigma2 <= 0.001) return(-100000)

  dat <- df_cluster

#   deliveryready_surface1 <- Smooth(df_spatial_hospital, sigma = param_sigma1)
#   deliveryready_surface2 <- Smooth(df_spatial_nonhospital, sigma = param_sigma2)

  get_surface_points <- function(surface, sigma_var, x_var, y_var) {
    surface1 <- Smooth(surface, sigma = sigma_var)
    get_pixel_value1 <- spatstat::as.function.im(surface1)
    # get pixel values (and fix where the function returns an empty vector)
    tmp <- mapply(get_pixel_value1, x = x_var, y = y_var)
    tmp[unlist(lapply(tmp, function(x) length(x)==0))] <- NA
    return(unlist(tmp))
  }

  dat$hosp <- get_surface_points(
    surface = df_spatial_hospital,
    sigma_var = param_sigma1,
    x_var = dat$xvar,
    y_var = dat$yvar )

  dat$nonhosp <- get_surface_points(
    surface = df_spatial_nonhospital,
    sigma_var = param_sigma2,
    x_var = dat$xvar,
    y_var = dat$yvar )

  dat <- dat %>%
    mutate(
      hosp = ifelse(is.na(hosp) | hosp < 0, 0, hosp),
      nonhosp = ifelse(is.na(nonhosp) | nonhosp < 0, 0, nonhosp),
      surface_var = mapply(max, hosp, nonhosp)
    )

  dat2 <- left_join(df_dhs, dat, by = c("clusterid" = "cluster_id")) %>%
    filter(complete.cases(.)) # figure out why 25 out of 1991 are NA

  y <- dat2$facilitybirth
  x <- cbind(dat2$wealth, dat2$educat, dat2$surface_var)

  os <- rep(1,length(x[,1]))
  x2 <- cbind(os,x)
  b <- param_regression[1:ncol(x2)]
  xb <- x2 %*% b
  sum(-1 * (y*log(1+exp(-xb)) + (1-y)*log(1+exp(xb))))
}


stval_names <- c(
  "constant", "wealth", "educat", "surfacevar", "sigma_hospital", "sigma_nonhospital" )

# coef_stval <- coef(glm(facilitybirth ~ wealth, data = df_dhs))
# stval <- c(coef_stval, 1, 1, 1) # constant, wealth, deliveryready, sigma1, sigma2
# stval <- c(-0.39, 0.126, 0.021, 0.65, 0.355) # dev variable
stval <- c(-0.39, 0.126, 2, 0.021, 40.0, 0.355)

result1 <- optim(
  par = stval, fn = llk.logit,
  method = "BFGS", hessian = T, control = list(fnscale = -1)
)

matrix(
  c(round(result1$par, digits = 4), round(sqrt(diag(solve(-1 * result1$hessian))), digits = 4)),
  dimnames = list(stval_names, c("Estimate", "Std. Error")),
  ncol = 2
)


tmp <- readRDS("data/results1.RDS")


# sig_val <- 2
sig_val <- result1$par[2]
anc_surface <- Smooth(df_spatial, sigma = sig_val)
plot(anc_surface, axes = TRUE)






