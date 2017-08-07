#
# 03_presentation_map.R
#

library(dplyr)
library(spatstat)

## maps for presentation

tmp9 <- readRDS("data/results_2sigmas_pca_and_all_others_v2.RDS")
result1 <- tmp9

# pixel_length <- 0.01 # takes longer; higher resolution
pixel_length <- 0.5

zoom <- TRUE

if (zoom) {

  xlo <- min(df$xvar) + 189 - pixel_length
  xhi <- min(df$xvar) + 203.9 + pixel_length
  ylo <- min(df$yvar) + 179 - pixel_length
  yhi <- min(df$yvar) + 194 + pixel_length

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

df_list3 <- lapply(df_list, function(x) {

  # x <- df_list[[2]] # dev variable

  tmp <- x %>% # keep points within (potentially zoomed) area
    filter(xvar >= xlo & xvar <= xhi) %>%
    filter(yvar >= ylo & yvar <= yhi) %>%
    dplyr::select(xvar, yvar, deliveryready)

  tmp2 <- rbind(tmp, dat_grid) # add grid of zeros

  tmp3 <- ppp( # create spatial data frame
    x = tmp2$xvar,
    y = tmp2$yvar,
    window = window1
  )

  marks(tmp3) <- tmp2$deliveryready # add facility values

  return(list(tmp3, max(tmp$deliveryready)))

})

(sigma_hosp <- result1$par[7])
(sigma_nonhosp <- result1$par[8])

s1 <- Smooth(df_list3[[1]][[1]], sigma = sigma_hosp)
s2 <- Smooth(df_list3[[2]][[1]], sigma = sigma_nonhosp)

s1$v <- s1$v * (df_list3[[1]][[2]] / max(s1$v))
s2$v <- s2$v * (df_list3[[2]][[2]] / max(s2$v))

plot(s1, axes = TRUE)
plot(s2, axes = TRUE)


tmp_v <- matrix(mapply(max, s1$v, s2$v),
  nrow = ncol(s1$v)
)

tmp_surface <- s1
tmp_surface$v <- tmp_v
tmp_surface$xrange <- c(0, xhi-xlo)
tmp_surface$yrange <- c(0, yhi-ylo)
tmp_surface$xcol <- tmp_surface$xcol - min(tmp_surface$xcol)
tmp_surface$yrow <- tmp_surface$yrow - min(tmp_surface$yrow)

plot(tmp_surface, main = " ", axes = TRUE)
mtext("Kilometers       ", side = 1, line = 2.4)


x_tmp <- seq(0, 2.5, by = 0.01)
y_hosp <- dnorm(x_tmp, mean = 0, sd = sigma_hosp)
y_nonhosp <- dnorm(x_tmp, mean = 0, sd = sigma_nonhosp)

plot(x_tmp, y_nonhosp, type = "n", xlab = "Kilometers", ylab = "Density", main = " ")
lines(x_tmp, y_nonhosp)
lines(x_tmp, y_hosp, lty = 2)
legend("topright",
  legend = c("Hospitals", "Clinics"),
  lty = c(2,1), cex = 0.9, bty = "n"
)

# x <- tmp_surface$xcol
# y <- tmp_surface$yrow
# val <- tmp_v
#
# library(raster)
#
# colors1 <- colorRampPalette(c("white", "darkblue"))
#
# tmp4 <- raster(ncols = length(x), nrows = length(y) )
# values(tmp4) <- as.vector(val)
# plot(tmp4)
# plot(tmp4, col = colors1(val))



