################################################################################
                  ##### Preparing environmental layers #####
################################################################################


# This script prepares environmental rasters for maxent models. Each raster is
# resampled to the same resolution and extent as the Great Britain land cover 
# 2000 map, which already matches the occurrence data for the hoverfly 
# recording scheme. 


# Load packages.
library(raster)
library(virtualspecies)




################################################################################
    #### Prepare UK Land cover data as the standard for other layers ####


# Load land use raster.
land_use <- raster("Data/Land Use/LCM2000_GB_1K_DOM_TAR.tif")

# Set coordinate reference system as UK OSGB 1936.
crs(land_use) <- '+init=EPSG:27700'

# Remove values that are considered sea. (This is different for newer land cover maps)
land_use[land_use==1] <- NA

# Change the name.
names(land_use) <- "Land_Use"

# Plot land use.
plot(land_use)




################################################################################
                        #### Add the climate data ####


# Create a list of file names to read into R of each climate layer.
bioclim_filenames <- list.files("Data/World Clim/WGS", full.names = T)

# Create a raster stack of climate variables.
bio_stack <- stack(lapply(bioclim_filenames, raster))

# Roughly crop to the extent of the UK in Long/Lat to make it quicker to reproject.
WGS_extent <- c(-8, 2, 49, 62)
bio_stack <- crop(bio_stack, WGS_extent)

# Project into OSGB coordinates.
bio_stack <- projectRaster(bio_stack, crs ='+init=EPSG:27700', method = 'bilinear')

# Resample rasters so they have the exact resolution as land use.
bio_stack <- resample(bio_stack, land_use, method="bilinear")

# Mask using land use.
bio_stack <- mask(bio_stack, land_use)

# Change the bioclim names to easier comprehend.
names(bio_stack) <- c("Annual_Mean_Temp",
                      "Mean_Diurnal_Range",
                      "Isothermality",
                      "Temp_Seasonality",
                      "Max_Temp_Warmest_Month",
                      "Min_Temp_Coldest_Month",
                      "Temp_Annual_Range",
                      "Mean_Temp_Wettest_Quarter",
                      "Mean_Temp_Driest_Quarter",
                      "Mean_Temp_Warmest_Quarter",
                      "Mean_Temp_Coldest_Quarter",
                      "Annual_Precip",
                      "Precip_Wettest_Month",
                      "Precip_Driest_Month",
                      "Precip Seasonality",
                      "Precip_Wettest_Quarter",
                      "Precip_Driest_Quarter",
                      "Precip_Warmest_Quarter",
                      "Precip_Coldest_Quarter",
                      "Elevation")

# Plot variables to check.
plot(bio_stack)

# Add land use to the stack.
bio_stack <- stack(land_use, bio_stack)

# Check for collinearity of variables. Use spearmans because 
# precipitation layers can be skewed. 
removeCollinearity(bio_stack, plot = TRUE,
                                   multicollinearity.cutoff = 0.7, 
                                   method="spearman")

# Create a directory for saving plots.
if (!dir.exists("Plots/Statistics")){
  dir.create("Plots/Statistics", recursive = TRUE)
}

# Use an edited function from virtual species that plots with spearmans label.
source("Code/functions.R")

# Plot the collinearity.
pdf("Plots/Statistics/collinearity.pdf")
remove_col_spearmans(bio_stack, plot = TRUE,
                   multicollinearity.cutoff = 0.7, 
                   method="spearman")
dev.off()

# Also save as a tiff.
tiff("Plots/Statistics/collinearity.tif")
remove_col_spearmans(bio_stack, plot = TRUE,
                     multicollinearity.cutoff = 0.7, 
                     method="spearman")
dev.off()

# Select the layers that aren't colinear. When layers are grouped, 
# we've selected layers based on biological relevance and simple interpretation,
# however because we model multiple empirical and virtual species with varied 
# niches, we do not interpret changes in variable importance during bias correction. 
env_stack <- subset(bio_stack, c("Land_Use",
                                 "Isothermality",
                                 "Temp_Seasonality",
                                 "Max_Temp_Warmest_Month",
                                 "Min_Temp_Coldest_Month",
                                 "Mean_Temp_Wettest_Quarter",
                                 "Mean_Temp_Driest_Quarter",
                                 "Annual_Precip",
                                 "Elevation"))


# Create a directory for saving layers.
if (!dir.exists("Data/Environmental Layers")){
  dir.create("Data/Environmental Layers", recursive = TRUE)
}

# Individually export rasters.
for (x in 1:9){
  filename <- paste("Data/Environmental Layers/", names(env_stack)[x], ".tif", sep="")
  writeRaster(env_stack[[x]], filename, overwrite=TRUE)
}





                       ##### End of Script. #####