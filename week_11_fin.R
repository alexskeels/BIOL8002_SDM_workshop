# If not installed use the pacman package to both install and load packages
install.packages("pacman")
library(pacman)

# Use pacman::p_load 
pacman::p_load(terra)  # For raster and vector data
pacman::p_load(sf)  # For vector data
pacman::p_load(sp) 
pacman::p_load(ggplot2)  # For plotting
pacman::p_load(rnaturalearth) # For polygons
pacman::p_load(rnaturalearthdata)
pacman::p_load(dplyr) # For data wrangling

# To make Alpha-Hulls we need to install a package from GitHub as well
remotes::install_github("babichmorrowc/hull2spatial")
pacman::p_load(hull2spatial)
pacman::p_load(alphahull)  # Required by hull2spatial internally

occ <- read.csv("distribution/frilled_lizard_ALA.csv")

# take a look at first few rows
head(occ)

# You can see there is a lot of meta-data - for example we can see what kinds of observations the records are from
unique(occ$basisOfRecord)

# how many records do we have?
nrow(occ)

# First, remove any rows with missing coordinates
occ <- occ %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude))

# Ssing the sf pacakge, turn the records in spatial points 
occ_sf <- st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Bounding box from your occurrence records - this is the extent the records cover
bbox <- st_bbox(occ_sf)

# Load the world map from rnaturalearth pacakge
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot using ggplot2
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray40") + # this line adds the world map
  geom_sf(data = occ_sf, aes(colour=basisOfRecord),  alpha = 0.6, size = 1) + # this line adds the occurence records
  coord_sf( # this section crops the map to the area surrounding the occurence records, ratehr than the whole globe
    xlim = c(bbox["xmin"]-10, bbox["xmax"]+10), # add 10 degrees on either side 
    ylim = c(bbox["ymin"]-10, bbox["ymax"]+10),
    expand = FALSE) +
  theme_minimal() + # this adds the theme elements to make it pretty
  labs(title = "Occurrence Records for the Frilled Lizard", x = "Longitude", y = "Latitude") # this adds the labels


# read in the IUCN polygon using the sf package
iucn_sf <- st_read("distribution/frilled_lizard.shp")


# plot with our occurence records
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray40") + # this line adds the world map
  geom_sf(data = occ_sf, aes(), color = "#0072B2", alpha = 0.6, size = 1) + # this line adds the occurence records
  geom_sf(data = iucn_sf, fill = "#009E73", color = "#007F5F", alpha = 0.4) +
  coord_sf( # this section crops the map to the area surrounding the occurence records, ratehr than the whole globe
    xlim = c(bbox["xmin"]-10, bbox["xmax"]+10), # add 10 degrees on either side 
    ylim = c(bbox["ymin"]-10, bbox["ymax"]+10),
    expand = FALSE) +
  theme_minimal() + # this adds the theme elements to make it pretty
  labs(title = "Occurrence Records for the Frilled Lizard", x = "Longitude", y = "Latitude") # this adds the labels



# choose buffer distance in meters
buffer_distance <- 100000 # CHANGE THIS

# create buffer in sf
buffer_iucn <- st_buffer(iucn_sf, buffer_distance)

# Buffers + original IUCN range + all occurrence points
ggplot() +
  geom_sf(data = buffer_iucn, color = "black", show.legend = TRUE, alpha=0.2) +
  geom_sf(data = iucn_sf, fill = "#009E73", color = "#007F5F", alpha = 0.5) +
  geom_sf(data = occ_sf, color = "#D55E00", alpha = 0.7, size = 1) +
  scale_fill_viridis_d(name = "Buffer size") +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(title = "IUCN Range with Buffer and Occurrence Points",
       x = "Longitude", y = "Latitude")


# subset only the points that fall within the polygon
occ_cropped <- occ_sf[lengths(st_within(occ_sf, buffer_iucn)) > 0, ]

# create a data frame and we'll save it for use in next week's prac
occ_cropped_df <- as.data.frame(st_coordinates(occ_cropped))

# write it as a csv
write.csv(occ_cropped_df, file="distribution/frilled_lizard_ALA_cropped.csv", row.names = F)

# crop the polygon by the coastline
sf_use_s2(FALSE)
buffer_iucn <- st_intersection(st_make_valid(buffer_iucn), world)
sf_use_s2(TRUE)

# Plot polygon and the cropped points
ggplot() +
  geom_sf(data = buffer_iucn, fill = "#56B4E9", color = "black", alpha = 0.3) +
  geom_sf(data = iucn_sf, fill = "#009E73", color = "#007F5F", alpha = 0.5) +
  geom_sf(data = occ_cropped, color = "#D55E00", size = 1) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(title = "Occurrence Points Within Buffered Polygon",
       x = "Longitude", y = "Latitude")


# Convert sf to data.frame to extract coordinates
coords <- st_coordinates(occ_cropped)
lon <- coords[, "X"]
lat <- coords[, "Y"]

# Compute Convex Hull indices
hull_indices <- chull(lon, lat)

# Ensure it's a closed polygon (repeat first point at the end)
hull_indices <- c(hull_indices, hull_indices[1])

# Create a matrix of hull coordinates
hull_coords <- cbind(lon[hull_indices], lat[hull_indices])

# Convert convex hull to an sf POLYGON
hull_polygon <- st_polygon(list(hull_coords))
hull_sf <- st_sfc(hull_polygon, crs = 4326)
hull_sf <- st_sf(geometry = hull_sf)

# crop the convex hull by the coastline
sf_use_s2(FALSE)
hull_sf <- st_intersection(st_make_valid(hull_sf), world)
sf_use_s2(TRUE)

# Plot Convex Hull
ggplot() +
  geom_sf(data = occ_cropped,  color = "#DD1C77", size = 0.8, alpha=0.5) +
  geom_sf(data = hull_sf, fill = NA, color = "#08519C") +
  
  theme_minimal() +
  labs(title = "Convex Hull from Cropped Occurrences")



alpha <- 5 # choose alpha value

# Compute Alpha Hull
xy_df <- data.frame(x = lon, y = lat)
xy_df <- unique(round(xy_df, 4))

ahull_obj <- ahull(xy_df, alpha = alpha)


# Convert alpha hull to sf POLYGON
ahull_spatial <- ahull2poly(ahull_obj)

# separate each sub-polygon
sub_polys <- slot(ahull_spatial, "polygons")

# Convert each sub-polygon to its own SpatialPolygons object
poly_list <- lapply(seq_along(sub_polys), function(i) {
  SpatialPolygons(list(sub_polys[[i]]), proj4string = CRS(proj4string(ahull_spatial)))
})

# Convert each to sf and combine into one sf object
sf_list <- lapply(poly_list, st_as_sf)
ahull_sf <- do.call(rbind, sf_list)

# Add an ID column (optional)
ahull_sf$id <- seq_len(nrow(ahull_sf))


st_crs(ahull_sf) <- 4326  # Set CRS

# crop the convex hull by the coastline
sf_use_s2(FALSE)
ahull_sf <- st_intersection(st_make_valid(ahull_sf), world)
sf_use_s2(TRUE)


ggplot() +
  geom_sf(data = occ_cropped,  color = "#DD1C77", size = 0.8, alpha=0.5) +
  geom_sf(data = ahull_sf, fill = NA, color = "#D94801") +
  theme_minimal() +
  coord_sf(expand = FALSE) +
  labs(
    title = "Frilled Lizard Distributions and Range Envelopes",
    subtitle = "IUCN range, buffer zone, convex and alpha hulls, and occurrence records",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    #legend.position = "none",
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.background = element_rect(fill = "white")
  )


st_write(ahull_sf, "distribution/frilled_lizard_alpha_hull.shp", append=F)


# Add a new column to each sf object for identifying it in the legend
ahull_sf$layer <- "Alpha Hull"
hull_sf$layer <- "Convex Hull"
iucn_sf$layer <- "IUCN Range"
buffer_iucn$layer <- "Buffer Zone"
occ_cropped$layer <- "Occurrences"

# Plot with color mapped to `layer`
ggplot() +
  geom_sf(data = occ_cropped, aes(color = layer), size = 0.8, alpha = 0.5) +
  geom_sf(data = ahull_sf, aes(color = layer), fill = NA) +
  geom_sf(data = buffer_iucn, aes(color = layer), fill = NA) +
  geom_sf(data = iucn_sf, aes(color = layer), fill = NA) +
  geom_sf(data = hull_sf, aes(color = layer), fill = NA) +
  scale_color_manual(
    name = "Layer",
    values = c(
      "Occurrences" = "#DD1C77",
      "Alpha Hull" = "#D94801",
      "Buffer Zone" = "#6BAED6",
      "IUCN Range" = "#238B45",
      "Convex Hull" = "#08519C"
    )
  ) +
  theme_minimal() +
  coord_sf(expand = FALSE) +
  labs(
    title = "Frilled Lizard Distributions and Range Envelopes",
    subtitle = "IUCN range, buffer zone, convex and alpha hulls, and occurrence records",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.background = element_rect(fill = "white")
  )

