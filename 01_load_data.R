#
# 01_load_data.R
#
# Reed Sorensen
# November 2016
#

require(foreign)
require(dplyr)

dir <- "C:/Users/rsoren/Documents/prog/repos/csss510_class_project/"

# DHS
df1 <- read.dbf(paste0(dir, "data/DHS GPS/HTGE61FL.dbf")) %>%
  select(
    id = DHSID,  cluster = DHSCLUST, year = DHSYEAR,
    latitude = LATNUM, longitude = LONGNUM
  )

# SPA
# SPAID = DHSCC&"SPA"&SPAYEAR&SPAFACID (with 8 digits) from survey documentation
#

df2 <- read.dbf(paste0(dir, "data/SPA GPS/HTGE6AFL.dbf")) %>%
  select(
    id = SPAID, facility_id = SPAFACID,
    latitude = LATNUM, longitude = LONGNUM
  )

# routine action in facility: blood pressure (cf106a)
# -- add a variable to the geospatial data,
#    so we have something to map

tmp_spa <- read.dta(paste0(dir,
  "data/spa_2013_haiti/Haiti 2013/htan6adtsr/HTAN6AFLSR.dta")) %>%
  select(facility_id = inv_id, bp_test = cf106a)

table(tmp_spa$facility_id %in% df2$facility_id)



# kernel density estimation

facility_dat <- read.csv(paste0(dir, "data/haitiancclean.csv"))

df3 <- left_join(df2, facility_dat, by = c("facility_id" = "facil")) %>%
  select(facility_id, latitude, longitude, ancscore)


# 'spatstat approach'
par(mfrow = c(1,1))

# remove NAs
df3 <- subset(df3, !is.na(ancscore))
df3 <- df3[c(3,4,14,24), ]
# df3 <- df3[sample(1:nrow(df3), size = 10), ]

x_pct_min <- 0
x_pct_max <- 1
y_pct_min <- 0
y_pct_max <- 1

x_range <- c(min(df3$longitude) - 0.05, max(df3$longitude) + 0.05)
y_range <- c(min(df3$latitude) - 0.05, max(df3$latitude) + 0.05)

(xmin <- (x_range[2] - x_range[1]) * x_pct_min + x_range[1])
(xmax <- (x_range[2] - x_range[1]) * x_pct_max + x_range[1])
(ymin <- (y_range[2] - y_range[1]) * y_pct_min + y_range[1])
(ymax <- (y_range[2] - y_range[1]) * y_pct_max + y_range[1])

# df3 <- df3 %>%
#   filter(longitude > xmin & longitude <= xmax) %>%
#   filter(latitude > ymin & latitude < ymax)


library(spatstat)
df4 <- ppp(
  x = df3$longitude,
  y = df3$latitude,
  xrange = c(min(df3$longitude) - 0.05, max(df3$longitude) + 0.05),
  yrange = c(min(df3$latitude) - 0.05, max(df3$latitude) + 0.05)
)

# plot(df4)
# plot(Kest(df4))
# plot(density(df4))

marks(df4) <- df3$ancscore
anc_sigma <- 0.20
anc_smooth <- Smooth(df4, sigma = anc_sigma)
anc_smooth$v[anc_smooth$v < 0] <- 0
anc_smooth$v[anc_smooth$v > 8] <- 8
plot(anc_smooth, main = "")
mtext(paste("ANC score; sigma =", anc_sigma), line = 0.1)

hmap <- readRDS("data/HTI_adm0.rds")
plot(hmap, add = TRUE)




















