# Load packages needed for spatial analysis
library(raster)
library(plyr)
library(dplyr)
library(dismo)
library(spatialEco)

# Set working directory to the sampling bias folder
setwd("Sampling Bias")

# Tell rJava where to find the java for developers folder
# (this is stored in a folder within the current directory)
Sys.setenv(JAVA_HOME="Java/jdk-13.0.1/")
library(rJava)




#### Prepare data for modelling ####

# Load climate file names
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

# Load all occurences
hoverflies <- read.csv("Data/Hoverfly Data/D_ALL_1983_2002.csv")
colnames(hoverflies)[1] <- "species"

# Total all species and filter for those over 1000 
sum_hoverflies <- plyr::count(hoverflies, "species")
over1000 <- filter(sum_hoverflies, freq > 1000)
print(over1000)

# Create an empty list for data frames of occurrences for each species
over1000_species <- list()

# Create a dataframe for each species that's over 1000
i = 1
for (hoverfly in over1000$species){
  name <- gsub(" ", "_", hoverfly)
  # Filter for specific species
  dataframe <- filter(hoverflies, species == hoverfly)
  assign(name, dataframe)
  # Save for future reference
  over1000_species[i] <- name
  i = i + 1
}




#### Read in bias files ####

# Change directory to the bias file folder
setwd("Bias Files/")

# Store bias file names and remove folders
bias_filenames <- list.files()
bias_filenames <- bias_filenames[-1]
bias_filenames <- bias_filenames[-3]

# Read in the bias files
for (x in bias_filenames){
  # Remove file extension
  name <- gsub(".tif", "", x)
  # Read raster and assign name
  raster <- raster(x)
  assign(name, raster)
}

# Load in distance restricted bias files
setwd("Buffers")
list.files()
for (x in list.files()){
  name <- gsub("_biasfile.tif", "", x)
  raster <- raster(x)
  assign(name, raster)
}

# Load in target group bias files
setwd("../Target Group/")
list.files()
for (x in list.files()){
  name <- gsub("_bias_file.tif", "", x)
  raster <- raster(x)
  assign(name, raster)
}




#### Create the first model that has sampling bias (no bias file) ####

# Change working directory back to main folder
setwd("../..")
getwd()

# Create a new directory for results
dir.create("Results")
dir.create("Results/Model Data")
dir.create("Results/Model Maps")

# Create directory for true bias data
dir.create("Results/Model Data/True Bias")
dir.create("Results/Model Maps/True Bias")

# Loop through each species and run the maxent models
for(x in 1:length(over1000_species)){
  # Get occurrences
  occurrences <- eval(parse(text=over1000_species[x]))[,2:3]
  # Create a file pathway
  file_pathway <- paste("Results/Model Data/True Bias/", over1000_species[x], sep="")
  dir.create(file_pathway)
  # Create the maxent model with 10,000 background points
  model <- maxent(climate_stack, occurrences, factors="landuse", path=file_pathway, nbg=10000, overwrite=T)
  # Save the model results
  data <- as.data.frame(model@results)
  write.csv(data, paste(file_pathway, "/model_results.csv", sep=""))
  # Create the predicted map of each model
  map <- predict(model, climate_stack)
  map <- raster.transformation(map, trans="norm")
  # Save the map
  name <- paste(over1000_species[x], "true_bias_pred", sep="_")
  file_pathway <- paste("Results/Model Maps/True Bias/", name, ".tif", sep="")
  writeRaster(map, file_pathway, overwrite=T)
}




#### Run the models that use the sampling effort bias file ####

# Make a function for extracting 10,000 background points from a bias file
background <- function(x) {
  bg_points <- xyFromCell(x, sample(which(!is.na(values(x))), 10000, prob = values(x)[!is.na(values(x))]))  
  return(bg_points)
}

# Create background points for each species using the bias file
background_points_list <- list()
for (x in 1:length(over1000_species)){
  background_points <- background(sample_effort_bias_file)
  name <- paste(over1000_species[x], "background", sep="_")
  assign(name, background_points)
  background_points_list[x] <- name
}

# Create directory for sample effort data
dir.create("Results/Model Data/Sample Effort")
dir.create("Results/Model Maps/Sample Effort")

# Run the maxent model
for(x in 1:length(over1000_species)){
  occurences <- eval(parse(text=over1000_species[x]))[,2:3]
  background_points <- eval(parse(text=background_points_list[x]))
  file_pathway <- paste("Results/Model Data/Sample Effort/", over1000_species[x], sep="")
  dir.create(file_pathway)
  model <- maxent(climate_stack, occurences, background_points, factors="landuse", path=file_pathway, overwrite=T)
  data <- as.data.frame(model@results)
  write.csv(data, paste(file_pathway, "/model_results.csv", sep=""))
  map <- predict(model, climate_stack)
  map <- raster.transformation(map, trans="norm")
  name <- paste(over1000_species[x], "sample_effort_pred", sep="_")
  file_pathway <- paste("Results/Model Maps/Sample Effort/", name, ".tif", sep="")
  writeRaster(map, file_pathway, overwrite=T)
}




#### Run the models that use the target group bias file ####


# Create background points for each species using the bias file
background_points_list <- list()
for (x in 1:length(over1000_species)){
  buffer <- eval(parse(text=paste(over1000_species[x], "_target_group", sep="")))
  background_points <- background(buffer)
  name <- paste(over1000_species[x], "background", sep="_")
  assign(name, background_points)
  background_points_list[x] <- name
}

# Create directory for target group data
dir.create("Results/Model Data/Target Group")
dir.create("Results/Model Maps/Target Group")

# Run the maxent model
for(x in 1:length(over1000_species)){
  occurences <- eval(parse(text=over1000_species[x]))[,2:3]
  background_points <- eval(parse(text=background_points_list[x]))
  file_pathway <- paste("Results/Model Data/Target Group/", over1000_species[x], sep="")
  dir.create(file_pathway)
  model <- maxent(climate_stack, occurences, background_points, factors="landuse", path=file_pathway, overwrite=T)
  data <- as.data.frame(model@results)
  write.csv(data, paste(file_pathway, "/model_results.csv", sep=""))
  map <- predict(model, climate_stack)
  map <- raster.transformation(map, trans="norm")
  name <- paste(over1000_species[x], "target_group_pred", sep="_")
  file_pathway <- paste("Results/Model Maps/Target Group/", name, ".tif", sep="")
  writeRaster(map, file_pathway, overwrite=T)
}




#### Run the models that use the buffer bias file ####

# Create background points for each species using the bias file
background_points_list <- list()
for (x in 1:length(over1000_species)){
  buffer <- eval(parse(text=paste(over1000_species[x], "_buffer", sep="")))
  background_points <- background(buffer)
  name <- paste(over1000_species[x], "background", sep="_")
  assign(name, background_points)
  background_points_list[x] <- name
}

# Create directory for buffer data
dir.create("Results/Model Data/Buffer")
dir.create("Results/Model Maps/Buffer")

# Run the maxent model
for(x in 1:length(over1000_species)){
  occurences <- eval(parse(text=over1000_species[x]))[,2:3]
  background_points <- eval(parse(text=background_points_list[x]))
  file_pathway <- paste("Results/Model Data/Buffer/", over1000_species[x], sep="")
  dir.create(file_pathway)
  model <- maxent(climate_stack, occurences, background_points, factors="landuse", path=file_pathway, overwrite=T)
  data <- as.data.frame(model@results)
  write.csv(data, paste(file_pathway, "/model_results.csv", sep=""))
  map <- predict(model, climate_stack)
  map <- raster.transformation(map, trans="norm")
  name <- paste(over1000_species[x], "Buffer_pred", sep="_")
  file_pathway <- paste("Results/Model Maps/Buffer/", name, ".tif", sep="")
  writeRaster(map, file_pathway, overwrite=T)
}




#### Run the models that use the population density bias file ####

# Create background points for each species using the bias file
background_points_list <- list()
for (x in 1:length(over1000_species)){
  background_points <- background(pop_dens_bias_file)
  name <- paste(over1000_species[x], "background", sep="_")
  assign(name, background_points)
  background_points_list[x] <- name
}

# Create directory for population density data
dir.create("Results/Model Data/Population Density")
dir.create("Results/Model Maps/Population Density")

# Run the maxent model
for(x in 1:length(over1000_species)){
  occurences <- eval(parse(text=over1000_species[x]))[,2:3]
  background_points <- eval(parse(text=background_points_list[x]))
  file_pathway <- paste("Results/Model Data/Population Density/", over1000_species[x], sep="")
  dir.create(file_pathway)
  model <- maxent(climate_stack, occurences, background_points, factors="landuse", path=file_pathway, overwrite=T)
  data <- as.data.frame(model@results)
  write.csv(data, paste(file_pathway, "/model_results.csv", sep=""))
  map <- predict(model, climate_stack)
  map <- raster.transformation(map, trans="norm")
  name <- paste(over1000_species[x], "pop_dens_pred", sep="_")
  file_pathway <- paste("Results/Model Maps/Population Density/", name, ".tif", sep="")
  writeRaster(map, file_pathway, overwrite=T)
}




#### Run the models that use the travel time bias file ####

# Create background points for each species using the bias file
background_points_list <- list()
for (x in 1:length(over1000_species)){
  background_points <- background(travel_time_bias_file)
  name <- paste(over1000_species[x], "background", sep="_")
  assign(name, background_points)
  background_points_list[x] <- name
}

# Create directory for travel time data
dir.create("Results/Model Data/Travel Time")
dir.create("Results/Model Maps/Travel Time")

# Run the maxent model
for(x in 1:length(over1000_species)){
  occurences <- eval(parse(text=over1000_species[x]))[,2:3]
  background_points <- eval(parse(text=background_points_list[x]))
  file_pathway <- paste("Results/Model Data/Travel Time/", over1000_species[x], sep="")
  dir.create(file_pathway)
  model <- maxent(climate_stack, occurences, background_points, factors="landuse", path=file_pathway, overwrite=T)
  data <- as.data.frame(model@results)
  write.csv(data, paste(file_pathway, "/model_results.csv", sep=""))
  map <- predict(model, climate_stack)
  map <- raster.transformation(map, trans="norm")
  name <- paste(over1000_species[x], "travel_pred", sep="_")
  file_pathway <- paste("Results/Model Maps/Travel Time/", name, ".tif", sep="")
  writeRaster(map, file_pathway, overwrite=T)
}

