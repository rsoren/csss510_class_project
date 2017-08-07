rm(list = ls())

library(geoR)
library(sp)
library(spatstat)
library(rgdal)
library(maptools)
library(splancs)
library(lattice)
library(SpatialEpi)


df <- read.csv("data/facility_data_anc.csv") %>%
  filter(!is.na(ancscore))

df[, c("xvar", "yvar")] <- latlong2grid(df[, c("longitude", "latitude")])

df <- df %>%
  mutate(
    xvar = xvar - min(xvar),
    yvar = yvar - min(yvar)
  )



zoom <- TRUE

if (zoom) {

  xlo <- 100
  xhi <- 120
  ylo <- 30
  yhi <- 50

} else {

  xlo <- min(df$xvar)
  xhi <- max(df$xvar)
  ylo <- min(df$yvar)
  yhi <- max(df$yvar)

}

pixel_length <- 0.5

x_points <- seq(xlo, xhi, by = pixel_length)
y_points <- seq(ylo, yhi, by = pixel_length)

df2 <- df %>%
  filter(xvar >= xlo & xvar <= xhi) %>%
  filter(yvar >= ylo & yvar <= yhi) %>%
  select(xvar, yvar, ancscore)

dat_grid <- expand.grid(xvar = x_points, yvar = y_points, ancscore = 0)

df3 <- rbind(df2, dat_grid)

window1 <- owin(
  xrange = c(min(df3$xvar), max(df3$xvar)),
  yrange = c(min(df3$yvar), max(df3$yvar)))

df_spatial <- ppp(
  x = df3$xvar,
  y = df3$yvar,
  window = window1
)

marks(df_spatial) <- df3$ancscore


sig_val <- 2
anc_surface <- Smooth(df_spatial, sigma = sig_val)
plot(anc_surface, axes = TRUE)

with(df2, points(xvar, yvar))















