# Load packages needed for spatial analysis
library(rgdal)
library(raster)
library(rgeos)
library(plyr)
library(dplyr)
library(ks)
library(spatialEco)

# Set working directory to the sampling bias folder
setwd("Sampling Bias")

# Load climate file names to mask against
climate_filenames <- paste("Data/World_clim/OSGB/", list.files("Data/World_clim/OSGB/"), sep="")

# Read in climate rasters
for (x in 1:length(climate_filenames)){
  name <- paste("world_clim", x, sep="")
  raster <- raster(climate_filenames[x])
  assign(name, raster)
}

# Load in the land use raster
land_use <- raster("Data/Land Use/landuse.tif")

# Stack cliamte rasters
climate_stack <- stack(land_use, world_clim1, world_clim2, world_clim3, world_clim4,
                       world_clim5, world_clim6, world_clim7, world_clim8, 
                       world_clim9, world_clim10, world_clim11, world_clim12, 
                       world_clim13, world_clim14, world_clim15, world_clim16, 
                       world_clim17, world_clim18, world_clim19)



#### Create true bias file based on sampling effort ####

# Load sampling effort data
sample_effort <- read.csv("Data/Sample_Effort/Sample_Effort_1983_2002.csv")

# Aggregate visits accross different years for each grid cell
sum_effort <- aggregate(sample_effort$nVisits, sample_effort[,3:4], FUN = sum)

# Create a matrix of coordinates
coords <- cbind(sum_effort[,1], sum_effort[,2])
head(coords)

# create scale so that sum of all sampling weights = total sample size
scale <- length(sum_effort[,3]) / sum(sum_effort[,3])
sum_effort[,4] <- sum_effort[,3] * scale

# Change names of columns
colnames(sum_effort) <- c("Northing", "Easting", "nVisits", "Scale")
head(sum_effort)

# Do a 2d kernel density estimate (using the scaled sampling effort as the weighting) 
# to create a probability surface of sampling effort
sample_density <- kde(coords, w=sum_effort$Scale)

# Create a raster of the sampling effort surface
sample_raster <- raster(sample_density)
plot(sample_raster)

# Define sampling effort raster as in OSGB
crs(sample_raster) <- '+init=EPSG:27700'

# Clip sample effort to same extent and resolution as land use data
sample_raster <- resample(sample_raster, land_use, method='bilinear')

# Mask sample effort to remove outside area
for (x in 1:20){
  sample_raster <- mask(sample_raster, climate_stack[[x]])
}

# Normalize bias file to between 0 and 1
sample_raster <- raster.transformation(sample_raster, trans="norm")

# Check that the plots match up
plot(land_use)
plot(sample_raster)
sample_raster

# Create a folder to store the bias files 
dir.create("Bias Files")

# Export the bias file
writeRaster(sample_raster, "Bias Files/sample_effort_bias_file.tif", overwrite=TRUE)




#### Create bias file based on target group data for each species ####
hoverflies <- read.csv("Data/Hoverfly Data/D_ALL_1983_2002.csv")
colnames(hoverflies)[1] <- "species"

# Total all species and filter for those over 1000 
sum_hoverflies <- plyr::count(hoverflies, "species")
over1000 <- filter(sum_hoverflies, freq > 1000)
print(over1000)

# Create a directory to save bias files
dir.create("Bias Files/Target Group")

# Loop through each species to create a unique bias file that's based on all other species in the HRS
for(hoverfly in over1000$species){
  # Remove target species
  records <- hoverflies %>% filter(hoverfly != species)
  # Aggregate by coordinates
  sum_records <- as.data.frame(dplyr::count(records, xx, yy))
  # Extract coordinates
  coords <- cbind(sum_records[,1], sum_records[,2])
  # Create a scale
  scale <- length(sum_records[,3]) / sum(sum_records[,3])
  sum_records[,4] <- sum_records[,3] * scale 
  # Do a 2d kernel density estimation
  target_density <- kde(coords, w=sum_records[,4])
  # Create raster
  target_raster <- raster(target_density)
  # Define in OSGB
  crs(target_raster) <- '+init=EPSG:27700'
  # Clip data to the same resolution/extent
  target_raster <- resample(target_raster, land_use, method='bilinear') 
  # Mask bias file
  for (x in 1:20){
    target_raster <- mask(target_raster, climate_stack[[x]])
  }
  # Normalize bias file between 0 and 1
  target_raster <- raster.transformation(target_raster, trans="norm")
  # Check each map
  plot(target_raster, main=hoverfly)
  # Create a filepathway and export the raster
  name <- gsub(" ", "_", hoverfly)
  pathway <- paste("Bias Files/Target Group/", name, "_target_group_bias_file.tif", sep="")
  writeRaster(target_raster, pathway, overwrite=TRUE)
}




#### Create buffer bias files ####

# Create a blank raster to use as a template for each bias file
blank_raster <- raster(ncol=700, nrow=1300, crs=crs(land_use))
extent(blank_raster) <- extent(land_use)

# Create a directory to store bias files
dir.create("Bias Files/Buffers")

# Loop through and create a bias file for each species with over 1000 records
for (hoverfly in over1000$species){
  # Create a species specific dataframe
  dataframe <- filter(hoverflies, species == hoverfly)
  # Create spatial points of all the occurrences
  points <- SpatialPoints(dataframe[,2:3], proj4string = CRS('+init=EPSG:27700'))
  # Rasterize the points using the blank raster
  raster <- rasterize(points, blank_raster, field=1)
  # Add a 5km buffer around each point
  raster_buffer <- buffer(raster, 10000)
  # Mask the raster
  for (x in 1:20){
    raster_buffer <- mask(raster_buffer, climate_stack[[x]])
  }
  # Plot each raster to check
  plot(raster_buffer, main=hoverfly)
  # Create a name for the pathway and export the bias file
  name <- gsub(" ", "_", hoverfly)
  pathway <- paste("Bias Files/Buffers/", name, "_buffer_biasfile.tif", sep="")
  writeRaster(raster_buffer, pathway, overwrite=TRUE)
}




#### Create a population density bias file ####

# Load in Gridded Human Population Density of the World
Pop_Dens_East_Europe <- raster("Data/Pop_dens/gpw_v4_population_density_rev11_2000_30_sec_3.asc")
plot(Pop_Dens_East_Europe)
Pop_Dens_West_Europe <- raster("Data/Pop_dens/gpw_v4_population_density_rev11_2000_30_sec_2.asc")
plot(Pop_Dens_West_Europe)

# Merge east and west europe into one raster
Europe <- raster::merge(Pop_Dens_West_Europe, Pop_Dens_East_Europe)

# Crop to roughly uk extent in WGS
WGS_extent <- c(-8, 2, 49, 62)
Pop_Dens_UK_WGS <- crop(Europe, WGS_extent)
plot(Pop_Dens_UK_WGS)

# Project into OSGB
Pop_Dens_UK <- projectRaster(Pop_Dens_UK_WGS, crs ='+init=EPSG:27700', method = 'bilinear')

# Clip map to UK resolution
Pop_Dens_UK <- resample(Pop_Dens_UK, land_use, method="bilinear")

# Mask map and plot
for (x in 1:20){
  Pop_Dens_UK <- mask(Pop_Dens_UK, climate_stack[[x]])
}
plot(Pop_Dens_UK)

# Values are logged to reduce overcompensation
UK_Log <- log(Pop_Dens_UK)

# '-inf' values from logging 0 values are set to NAs
remove_infs <- function(x) { x[x==-Inf] <- NA; return(x)}
UK_Log <- calc(UK_Log, remove_infs)
summary(UK_Log)

#Values are shifted to a positive number (min value is approx -3.1)
UK_Log_Pos <- UK_Log - minValue(UK_Log)

#normalize bias file to between 0 and 1
UK_Log_Pos <- raster.transformation(UK_Log_Pos, trans="norm")
plot(UK_Log_Pos)

#export bias file
writeRaster(UK_Log_Pos, "Bias Files/pop_dens_bias_file.tif", overwrite=TRUE)




#### Create a bias file from travel time ####

# Load in Gridded Travel Time of the World
Travel_Time <- raster("Data/Travel_time/acc_50k.tif")
plot(Travel_Time)

# Crop to roughly the UK
Travel_Time_UK <<- crop(Travel_Time, WGS_extent)

# Project into OSGB
Travel_Time_UK <- projectRaster(Travel_Time_UK, crs ='+init=EPSG:27700', method = 'bilinear')
plot(Travel_Time_UK)

# Clip map to UK resolution
Travel_Time_UK <- resample(Travel_Time_UK, land_use, method="bilinear")

# Mask map
for (x in 1:20){
  Travel_Time_UK <- mask(Travel_Time_UK, climate_stack[[x]])
}
plot(Travel_Time_UK)

# Values are logged to reduce overcompensation
Travel_Time_Log <- log(Travel_Time_UK)

# '-inf' values from logging 0 values are set to NAs
Travel_Time_Log <- calc(Travel_Time_Log, remove_infs)
summary(Travel_Time_Log)
plot(Travel_Time_Log)

# values are negated to revserse the relationship (so that the most remove locations have the lowest probability of being selected)
Travel_Time_Log <- Travel_Time_Log * -1

# Values are shifted to a positive number (min value is approx -6.6)
Travel_Time_Log_Pos <- Travel_Time_Log - minValue(Travel_Time_Log)

# Normalize bias file to between 0 and 1
Travel_Time_Log_Pos <- raster.transformation(Travel_Time_Log_Pos, trans="norm")
plot(Travel_Time_Log_Pos)

# Export bias file
writeRaster(Travel_Time_Log_Pos, "Bias Files/travel_time_bias_file.tif", overwrite=TRUE)

