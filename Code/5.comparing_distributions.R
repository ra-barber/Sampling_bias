# Load packages needed for spatial analysis
library(raster)
library(spatialEco)
library(stringr)
library(dplyr)
library(ENMTools)

# Set working directory to the sampling bias folder
setwd("Sampling Bias")

#### Set up the dataframe and rasters ####

# Make a dataframe for storing results
results <- data.frame(method = character(), 
                      species = character(), 
                      rastername=character(), 
                      samplenumber=numeric(),
                      Schoeners = numeric(),
                      Hellingers = numeric(),
                      Spearmans = numeric(),
                      centroidX = numeric(),
                      centroidY = numeric(),
                      changeinX = numeric(),
                      changeinY = numeric(),
                      centroidshift = numeric(),
                      rangesize = numeric(),
                      stringsAsFactors=FALSE)


# Generate the filepathways of all the maps
map_filepathways <- list.files("Results/Model Maps", recursive = T)

# Read in the rasters and edit their names to alphabetically sort in a logical order
for (x in 1:length(map_filepathways)){
  # Read in each raster
  raster <- raster(paste("Results/Model Maps", map_filepathways[x], sep="/"))
  # Strip the file extensions
  name <- gsub(".tif", "", map_filepathways[x])
  # Replaces underscores after the species name with a letter so plots are in order when checking
  name <- str_replace(name, "_(?=sample_effort_)", "_A_")
  name <- str_replace(name, "_(?=true)", "_B_")
  name <- str_replace(name, "_(?=target)", "_C_")
  name <- str_replace(name, "_(?=Buffer)", "_D_")
  name <- str_replace(name, "_(?=pop)", "_E_")
  name <- str_replace(name, "_(?=travel)", "_F_")
  # Remove the parent folders from names
  name <- str_remove(name, ".*/")
  # Assign rasters with names
  assign(name, raster)
  
  # Strip extra information from raster names
  species <- str_remove(name, "_[A-F]_.*")
  method <- str_remove(name, ".*_[A-F]_")
  method <- str_remove(method, "_pred")
  
  # Store information on species, method of bias correction, and the name of the raster
  results[x,1] <- method
  results[x,2] <- species
  results[x,3] <- name
}

# Rearrange rasters by species and name
results <- arrange(results, species, rastername)




#### Compare distributions using metrics for range size, centroid position, and niche overlap ####

# Create a list of species
species_list <- base::unique(results$species)

# Define a matrix that reclassify uses to turn continous maps into binary maps with a threshold of 0.5
classes <- matrix(c(0,0.5,0,0.5,1,1), nrow=2,byrow=T)

# Loop through each species, and then through each bias correction method, to generate metrics
for (x in species_list){
  # Filter results for a single species
  single_species <- filter(results, species==x)
  # Loop through each method
    for (y in 1:6){
      # Pull out the specific raster using the raster name
      raster_name <- single_species[y,3]
      raster <- eval(parse(text=raster_name))
      # Get the sample effort distribution for the species
      sample_effort <- eval(parse(text=single_species[1,3]))
      # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
      overlaps <- raster.overlap(sample_effort, raster)
      # Turn the raster into points for calculating the centroid (of specfic method and sample effort for the species)
      raster_points <- rasterToPoints(raster, spatial = T)
      sample_points <- rasterToPoints(sample_effort, spatial = T)
      # Calcuate raster centroid and save as a spatial point (Package: spatialEco)
      raster_centroid <- wt.centroid(raster_points)
      sample_centroid <- wt.centroid(sample_points)
      # Record the position of the centroid and sample effort centroid
      centroidX <- xmin(raster_centroid)
      centroidY <- ymin(raster_centroid)
      sampleX <- xmin(sample_centroid)
      sampleY <- ymin(sample_centroid)
      # Calculate the change in Northing and Easting
      change_in_x <- centroidX - sampleX
      change_in_y <- centroidY - sampleY
      # Calculate the Euclidean distance
      magnitude <- sqrt((change_in_x^2)+(change_in_y^2))
      # Use the predefined matrix to turn rasters in binary maps with 0.5 threshold
      raster_bin <- reclassify(raster, classes)
      sample_bin <- reclassify(sample_effort, classes)
      # Calculate the number of cells, and thus area of cells with over 0.5 habitat suitability
      raster_area <- ncell(raster_bin[raster_bin==1])
      sample_area <- ncell(sample_bin[sample_bin==1])
      # Get the relative change in area to account for small and large distributions
      range_size <- raster_area / sample_area
      # Converts so that increases are centred at zero
      if (range_size > 1){
        range_size = range_size - 1
      } else {
        range_size = (1-range_size)*-1
      }
      # Convert scale to 100 for percentages
      range_size = range_size * 100
      # Add the metrics to results
      results[results$rastername==raster_name,5] <- overlaps$D
      results[results$rastername==raster_name,6] <- overlaps$I
      results[results$rastername==raster_name,7] <- overlaps$rank.cor
      results[results$rastername==raster_name,8] <- centroidX
      results[results$rastername==raster_name,9] <- centroidY
      results[results$rastername==raster_name,10] <- change_in_x
      results[results$rastername==raster_name,11] <- change_in_y
      results[results$rastername==raster_name,12] <- magnitude
      results[results$rastername==raster_name,13] <- range_size
    }
}




#### A loop to add in the sample number for each method (will be the same if there are no bugs) ####

# Get a list of method types as recorded in directories
method_types <- list.files("Results/Model Data")

# Change the method names to match directories
results$method <- c("Sample Effort", "True Bias", "Target Group", "Buffer", "Population Density", "Travel Time")

# Remove the period from species names to match directories
results$species <- str_remove(results$species, "\\.")

# Loop through the method types and add each sample number to results
for (x in method_types){
  # Get the pathway for each bias correction method
  method_pathway <- paste("Results/Model Data", x, sep="/")
  # Get the species list as recorded in directories
  species_list <- list.files(method_pathway)
  # Loop through species in each bias correction method
  for (y in species_list){
    # Get the pathway for the results table
    species_pathway <- paste(method_pathway, y, "model_results.csv", sep="/")
    # Read in the csv
    table <- read.csv(species_pathway)
    # Save the sample number in the row with matching method and species
    results$samplenumber[results$method==x & results$species==y] <- table[1,2]
  }
}





# Export results
dir.create("Results/Statistics")
write.csv(results, "Results/Statistics/comparison_metrics.csv", row.names = F)


