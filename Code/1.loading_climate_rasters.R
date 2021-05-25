# Load packages needed for spatial analysis
library(rgdal)
library(raster)

# Set working directory to the sampling bias folder
setwd("Sampling Bias")




#### Prepare UK Land Use data which will be used as the standard dimensions for all other layers ####

# Load Landuse raster
land_use <- raster("Data/Land Use/LCM2000_GB_1K_DOM_TAR.tif")

# Set coordinate reference system as UK OSGB 1936
crs(land_use) <- '+init=EPSG:27700'

# Check axis in plot to check it's worked okay
plot(land_use)

# Remove values that are considered sea
land_use[land_use==1] <- NA
summary(land_use)

# Check that the sea has been removed
plot(land_use)




#### Add the climate data ####

# Create a list of filenames to read into R of each climate layer
bioclim_filenames <- list.files("Data/World_clim/WGS/")

# Create a blank list to store the climate raster names
clim_rasters <- list()

#load in the rasters by looping through the filenames
i = 1
for (file_name in bioclim_filenames){
  #add the directory pathway
  filepathway <- paste("Data/World_clim/WGS/", file_name, sep = "")
  #read in the raster
  raster <- raster(filepathway)
  #create a name for each one 
  name <- paste("BioClim", i, sep="")
  #save the name in a list for easy referencing
  clim_rasters[i] <- name
  #assign the name to the raster
  assign(name, raster)
  #add the counter
  i = i + 1
}




#### Crop world climate maps to roughly the UK ####

# Roughly define the extent of the UK in Long/Lat to make it easier to project coordinate systems
WGS_extent <- c(-8, 2, 49, 62)

# Remove excess data from climate rasters
for(i in 1:19){
  # Get each raster name
  name <- clim_rasters[[i]]
  # Pull out the raster object
  raster <- eval(parse(text=name))
  # Crop the raster to roughly the UK extent
  raster <- crop(raster, WGS_extent)
  # Assign each raster back to the name
  assign(name, raster)
}

# Plot maps to check they were all cropped properly
for (i in 1:19){
  name <- clim_rasters[[i]]
  raster <- eval(parse(text=name))
  plot(raster, main=name)
}




#### Project Maps into OSGB Coordinates and match extent/resolution to Land Use ####

# Create a new list for OSGB maps
OSGB_names <- list()

# Project into the new coordinate system
for(i in 1:19){
  name <- clim_rasters[[i]]
  raster <- eval(parse(text=name))
  # Project each raster into OSGB 1936 using the bilinear method
  raster <- projectRaster(raster, crs ='+init=EPSG:27700', method = 'bilinear')
  # Create new names
  new_name <- paste(name, "_OSGB", sep="")
  # Assign raster
  assign(new_name, raster)
  # Save name to list for easy referencing
  OSGB_names[i] <- new_name
}


# Plot new maps to check they were projected to OSGB properly
for (i in 1:19){
  name <- OSGB_names[[i]]
  raster <- eval(parse(text=name))
  plot(raster, main=name)
}

# Resample OSGB maps to OSGB 1KM cell resolution to match hoverfly data resolution
for (i in 1:19){
  name <- OSGB_names[[i]]
  raster <- eval(parse(text=name))
  # Use the Land Use map as a reference for the correct dimensions
  raster <- resample(raster, land_use, method='bilinear')
  assign(name, raster)
}

# Mask climate maps so that they match the Land Use map (remove Ireland/France)
for (i in 1:19){
  name <- OSGB_names[[i]]
  raster <- eval(parse(text=name))
  # Mask sets all excess cells to NA
  raster <- mask(raster, land_use)
  assign(name, raster)
}

# Check that maps have come out okay
for (i in 1:length(clim_rasters)){
  name <- OSGB_names[[i]]
  raster <- eval(parse(text=name))
  plot(raster, main=name)
}




#### Export Maps ####

# Create directoy for OSGB maps
dir.create("Data/World_clim/OSGB/")

# Export Maps
for (i in 1:length(clim_rasters)){
  name <- OSGB_names[[i]]
  raster <- eval(parse(text=name))
  # Create file pathway to save to
  filepathway <- paste("Data/World_clim/OSGB/", name, ".tif", sep="")
  # Write Raster
  writeRaster(raster, filepathway, overwrite=TRUE)
}

writeRaster(land_use, "Data/Land Use/landuse.tif")

