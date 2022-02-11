################################################################################
                     ##### Creating maxent models #####
################################################################################


# This script performs all the maxent species distribution models. First we 
# create two reference set of maps to compare different bias correction 
# techniques against. For empirical species, these are the maps generated using
# the sample effort bias file. For virtual species, these are the maps generated
# using maxent models with a uniform bias file and unbiased occurrences. After
# generating maps for the main analysis, we then generate similarity metrics 
# for virtual species against their known probability distributions which are 
# used to generate similarity metrics shown in the supplementary materials.


# Load packages.
library(raster)
library(dplyr)
library(spatialEco)
library(virtualspecies)
library(foreach)
library(doParallel)
# Set the pathway for Java.
Sys.setenv(JAVA_HOME="Java/jdk-13.0.1/")
library(rJava)
library(dismo)
library(ENMTools)
library(biomod2)





################################################################################
                        ##### Read in Data #######


# Read in environmental rasters.
bioclim_filenames <- list.files("Data/Environmental Layers/", full.names = T)
bio_stack <- stack(lapply(bioclim_filenames, raster))

# Read in over 1000 occurrence data.
hoverflies <- read.csv("Data/Hoverfly Data/D_ALL_1983_2002.csv")
hoverfly_species <- unique(hoverflies$Species)

# Read in virtual occurrences.
virtual_occurences <- read.csv("Data/Virtual/Occurrences/all_occurrences.csv")
virtual_species <- unique(virtual_occurences$random_name)

# Read in the uniform bias file.
uniform_biasfile <- raster("Data/Bias Files/Uniform_biasfile.tif")

# Read in the sample effort bias file.
sample_effort_biasfile <- raster("Data/Bias Files/SampleEffort_biasfile.tif")

# Read in the population density bias file.
pop_dens_biasfile <- raster("Data/Bias Files/PopDens_biasfile.tif")

# Read in the sample effort bias file.
travel_time_biasfile <- raster("Data/Bias Files/TravelTime_biasfile.tif")

# Read in the target group bias files.
target_biasfile <- raster("Data/Bias Files/TargetGroup_biasfile.tif")

# Read in the buffer 10k bias files.
buffer_filenames <- list.files("Data/Bias Files/Buffer 10k/", full.names = T)
buffer10k_stack <- stack(lapply(buffer_filenames, raster))

# Read in the virtual buffer 10k bias files.
buffer_filenames <- list.files("Data/Bias Files/Virtual Buffer 10k/", full.names = T)
virtual_buffer10k_stack <- stack(lapply(buffer_filenames, raster))

# Set the seed.
set.seed(12345)




################################################################################
                     ##### Create reference models #####

tictoc::tic()
# Set up parallel computing.
registerDoParallel(cores = 8)

# Loop through virtual species and create the models using unbiased presence-only
# occurrence data.
virtual_reference <- foreach(x=1:length(virtual_species),
                  .inorder = TRUE,
                  .errorhandling = "stop",
                  .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- virtual_species[x]
    
    # Get the model name.
    model_name <- paste(species, "Unbiased", sep="_")
    
    # Filter occurrences for species.
    occurrences <- virtual_occurences %>% filter(random_name == species) %>% select(PO_x_unbiased, PO_y_unbiased)
  
    # Select 10000 background points without using probability.
    background_points <- randomPoints(mask = uniform_biasfile, n = 10000, prob = TRUE)
  
    # Create the file pathway.
    file_pathway <- paste("Results/Virtual/Unbiased",  species, sep="/")
  
    # Create the maxent model with 10,000 background points
    model <- maxent(bio_stack, occurrences, background_points, factors="Land_Use", 
                  path=file_pathway, nbg=10000, overwrite=TRUE, args=c("maximumiterations=2000"))
    
    # Save the model as an RDS object.
    save(model, file = paste(file_pathway, "/", model_name, "_model.rds", sep=""))
  
    # Create the predicted map of each model.
    map <- predict(model, bio_stack, args=c('outputformat=cloglog'))
    #map <- raster.transformation(map, trans = "norm")
  
    # Save the map.
    writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  }


# Loop through real species and create the models using sample effort bias file.
emperical_reference <- foreach(x=1:length(hoverfly_species),
                            .inorder = TRUE,
                            .errorhandling = "stop",
                            .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- hoverfly_species[x]
    
    # Get the model name.
    model_name <- paste(species, "SampleEffort", sep="_")
  
    # Filter occurrences for species.  #### Change selection.
    occurrences <- hoverflies %>% filter(Species == species) %>% select(Easting, Northing)
  
    # Select 10000 background points without using probability.
    background_points <- randomPoints(mask = sample_effort_biasfile, n = 10000, prob = TRUE)
  
    # Create the file pathway.
    file_pathway <- paste("Results/Emperical/SampleEffort",  species, sep="/")
  
    # Create the maxent model with 10,000 background points.
    model <- maxent(bio_stack, occurrences, background_points, factors="Land_Use", 
                  path=file_pathway, nbg=10000, overwrite=TRUE, args=c("maximumiterations=2000"))
  
    # Save the model as an RDS object.
    save(model, file = paste(file_pathway, "/", model_name, "_model.rds", sep=""))
  
    # Create the predicted map of each model.
    map <- predict(model, bio_stack, args=c('outputformat=cloglog'))
    #map <- raster.transformation(map, trans = "norm")
  
    # Save the map.
    writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  
}


tictoc::toc()

# Turn the lists into stacks.
virtual_stack <- stack(virtual_reference)
emperical_stack <- stack(emperical_reference)




################################################################################
                    ##### Model Virtual Species #####


# Load in function for modeling virtual species and recording metrics.
source("Code/functions.R")

# Model virtual species and extract metrics.
tictoc::tic()

# For each method, use the virtual maxent function from the functions.R script to
# create a maxent model and extract similarity metrics compared to models made
# with the unbiased virtual occurrences.
virtual <- foreach(x=1:length(virtual_species),
                  .inorder = FALSE,
                  .errorhandling = "stop",
                  .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    species <- virtual_species[x]
    occurrences <- virtual_occurences %>% filter(random_name == species) %>% select(PO_x_biased, PO_y_biased)
    virtual_maxent(species, "NoCorrection", uniform_biasfile)
    virtual_maxent(species, "SampleEffort", sample_effort_biasfile)
    virtual_maxent(species, "TargetGroup", target_biasfile)
    
    buffer_raster <- virtual_buffer10k_stack[[paste(species, "_Buffer10k_biasfile", sep="")]]
    virtual_maxent(species, "Buffer10k", buffer_raster)
    
    virtual_maxent(species, "PopDens", pop_dens_biasfile)
    virtual_maxent(species, "TravelTime", travel_time_biasfile)
  }

tictoc::toc()


tictoc::tic()
# Repeat models using emperical data.
emperical <- foreach(x=1:length(hoverfly_species),
                   .inorder = FALSE,
                   .errorhandling = "stop",
                   .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    species <- hoverfly_species[x]
    
    emperical_maxent(species, "NoCorrection", uniform_biasfile)
    emperical_maxent(species, "TargetGroup", target_biasfile)
    
    buffer_raster <- buffer10k_stack[[paste(species, "_Buffer10k_biasfile", sep="")]]
    emperical_maxent(species, "Buffer10k", buffer_raster)
    
    emperical_maxent(species, "PopDens", pop_dens_biasfile)
    emperical_maxent(species, "TravelTime", travel_time_biasfile)
  }

tictoc::toc()




################################################################################
     #### Create virtual distributions used in supplementary materials ####


# Read in probability of occurrences.
prob_filenames <- list.files("Data/Virtual/Probabilites/", full.names = T)
probability_stack <- stack(lapply(prob_filenames, raster))

# Create a biased set of virtual models.
biased_reference <- foreach(x=1:length(virtual_species),
                             .inorder = TRUE,
                             .errorhandling = "stop",
                             .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- virtual_species[x]
    
    # Get the model name.
    model_name <- paste(species, "Biased", sep="_")
    
    # Filter occurrences for species.
    occurrences <- virtual_occurences %>% filter(random_name == species) %>% select(PO_x_biased, PO_y_biased)
    
    # Select 10000 background points without using probability.
    background_points <- randomPoints(mask = uniform_biasfile, n = 10000, prob = TRUE)
    
    # Create the file pathway.
    file_pathway <- paste("Results/Virtual/Biased",  species, sep="/")
    
    # Create the maxent model with 10,000 background points
    model <- maxent(bio_stack, occurrences, background_points, factors="Land_Use", 
                    path=file_pathway, nbg=10000, overwrite=TRUE, args=c("maximumiterations=2000"))
    
    # Save the model as an RDS object.
    save(model, file = paste(file_pathway, "/", model_name, "_model.rds", sep=""))
    
    # Create the predicted map of each model.
    map <- predict(model, bio_stack, args=c('outputformat=cloglog'))
    #map <- raster.transformation(map, trans = "norm")
    
    # Save the map.
    writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  }


# Create a set of models using the biased data and the sample effort bias file.
biasfile_reference <- foreach(x=1:length(virtual_species),
                            .inorder = TRUE,
                            .errorhandling = "stop",
                            .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- virtual_species[x]
    
    # Get the model name.
    model_name <- paste(species, "BiasFile", sep="_")
    
    # Filter occurrences for species.
    occurrences <- virtual_occurences %>% filter(random_name == species) %>% select(PO_x_biased, PO_y_biased)
    
    # Select 10000 background points without using probability.
    background_points <- randomPoints(mask = sample_effort_biasfile, n = 10000, prob = TRUE)
    
    # Create the file pathway.
    file_pathway <- paste("Results/Virtual/BiasFile",  species, sep="/")
    
    # Create the maxent model with 10,000 background points
    model <- maxent(bio_stack, occurrences, background_points, factors="Land_Use", 
                    path=file_pathway, nbg=10000, overwrite=TRUE, args=c("maximumiterations=2000"))
    
    # Save the model as an RDS object.
    save(model, file = paste(file_pathway, "/", model_name, "_model.rds", sep=""))
    
    # Create the predicted map of each model.
    map <- predict(model, bio_stack, args=c('outputformat=cloglog'))
    #map <- raster.transformation(map, trans = "norm")
    
    # Save the map.
    writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  }


# Create raster stacks for both set of models.
biased_stack <- stack(biased_reference)
biasfile_stack <- stack(biasfile_reference)




################################################################################
    #### Create supplementary metrics to validate sampling effort ####


# Extract similarity metrics comparing the virtual species with unbiased data 
# to the probability distribution of the virtual species.
unbiased <- foreach(x=1:length(virtual_species),
                     .inorder = FALSE,
                     .errorhandling = "stop",
                     .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- virtual_species[x]
    
    # Create the model name.
    model_name <- paste(species, "Unbiased", sep="_")
    
    # Create the file pathway.
    file_pathway <- paste("Results/Virtual/Unbiased", species, sep="/")
    
    # Pull out the distributions we're testing.
    map <- virtual_stack[[paste(species, "_Unbiased_map", sep="")]]
    prob_raster <- probability_stack[[paste("prob_", species, sep="")]]
    
    # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
    prob_overlaps <- raster.overlap(prob_raster, map)
    
    # Calculate the centroids.
    prob_centroid <- wt.centroid(rasterToPoints(prob_raster, spatial = T))
    test_centroid <- wt.centroid(rasterToPoints(map, spatial = T))
    
    # Calculate differences.
    prob_change_in_x <- xmin(test_centroid) - xmin(prob_centroid)
    prob_change_in_y <- ymin(test_centroid) - ymin(prob_centroid)
    
    # Calculate the Euclidean distance in km.
    prob_magnitude <- sqrt((prob_change_in_x^2)+(prob_change_in_y^2))/1000
    
    # Calculate the area based on a binary classification.
    test_bin <- BinaryTransformation(map, 0.5)
    prob_bin <- BinaryTransformation(prob_raster, 0.5)
    
    # Calculate range size change.
    prob_change <- as.data.frame(BIOMOD_RangeSize(prob_bin, test_bin)$Compt.By.Models)
    
    # Add the metrics to results.
    similarity_metrics <- data.frame(method = "Unbiased", 
                                 species = species,
                                 
                                 # Niche Similarity.
                                 p_schoeners = prob_overlaps$D,
                                 p_hellingers = prob_overlaps$I,
                                 p_spearmans = prob_overlaps$rank.cor,
                                 
                                 # Centroid shift.
                                 p_centroidshift = prob_magnitude,
                                 
                                 # Range size changes.
                                 p_loss = prob_change$PercLoss,
                                 p_gain = prob_change$PercGain,
                                 p_total_change = prob_change$SpeciesRangeChange,
                                 p_range_size = prob_change$CurrentRangeSize,

                                 stringsAsFactors=FALSE)

write.csv(similarity_metrics, paste(file_pathway, "/", model_name, "_metrics.csv", sep=""), row.names = FALSE)
}


# Do the same for biased.
biased <- foreach(x=1:length(virtual_species),
                  .inorder = FALSE,
                  .errorhandling = "stop",
                  .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- virtual_species[x]
    
    # Create the model name.
    model_name <- paste(species, "Biased", sep="_")
    
    # Create the file pathway.
    file_pathway <- paste("Results/Virtual/Biased", species, sep="/")
    
    # Pull out the distributions we're testing.
    map <- biased_stack[[paste(species, "_Biased_map", sep="")]] # Using biased stack.
    prob_raster <- probability_stack[[paste("prob_", species, sep="")]]
    
    # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
    prob_overlaps <- raster.overlap(prob_raster, map)
    
    # Calculate the centroids.
    prob_centroid <- wt.centroid(rasterToPoints(prob_raster, spatial = T))
    test_centroid <- wt.centroid(rasterToPoints(map, spatial = T))
    
    # Calculate differences.
    prob_change_in_x <- xmin(test_centroid) - xmin(prob_centroid)
    prob_change_in_y <- ymin(test_centroid) - ymin(prob_centroid)
    
    # Calculate the Euclidean distance in km.
    prob_magnitude <- sqrt((prob_change_in_x^2)+(prob_change_in_y^2))/1000
    
    # Calculate the area based on a binary classification.
    test_bin <- BinaryTransformation(map, 0.5)
    prob_bin <- BinaryTransformation(prob_raster, 0.5)
    
    # Calculate range size change.
    prob_change <- as.data.frame(BIOMOD_RangeSize(prob_bin, test_bin)$Compt.By.Models)
    
    # Add the metrics to results.
    similarity_metrics <- data.frame(method = "Biased", 
                                     species = species,
                                     
                                     # Niche Similarity.
                                     p_schoeners = prob_overlaps$D,
                                     p_hellingers = prob_overlaps$I,
                                     p_spearmans = prob_overlaps$rank.cor,
                                     
                                     # Centroid shift.
                                     p_centroidshift = prob_magnitude,
                                     
                                     # Range size changes.
                                     p_loss = prob_change$PercLoss,
                                     p_gain = prob_change$PercGain,
                                     p_total_change = prob_change$SpeciesRangeChange,
                                     p_range_size = prob_change$CurrentRangeSize,
                                     
                                     stringsAsFactors=FALSE)
    
    write.csv(similarity_metrics, paste(file_pathway, "/", model_name, "_metrics.csv", sep=""), row.names = FALSE)
  }


# And the same for the bias file stack.
biasfile <- foreach(x=1:length(virtual_species),
                  .inorder = FALSE,
                  .errorhandling = "stop",
                  .packages = c("dismo", "raster", "spatialEco", "rJava", "dplyr", "ENMTools", "biomod2")) %dopar% 
  {
    # Pull out the species.
    species <- virtual_species[x]
    
    # Create the model name.
    model_name <- paste(species, "BiasFile", sep="_")
    
    # Create the file pathway.
    file_pathway <- paste("Results/Virtual/BiasFile", species, sep="/")
    
    # Pull out the distributions we're testing.
    map <- biasfile_stack[[paste(species, "_BiasFile_map", sep="")]] # Using biased stack.
    prob_raster <- probability_stack[[paste("prob_", species, sep="")]]
    
    # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
    prob_overlaps <- raster.overlap(prob_raster, map)
    
    # Calculate the centroids.
    prob_centroid <- wt.centroid(rasterToPoints(prob_raster, spatial = T))
    test_centroid <- wt.centroid(rasterToPoints(map, spatial = T))
    
    # Calculate differences.
    prob_change_in_x <- xmin(test_centroid) - xmin(prob_centroid)
    prob_change_in_y <- ymin(test_centroid) - ymin(prob_centroid)
    
    # Calculate the Euclidean distance in km.
    prob_magnitude <- sqrt((prob_change_in_x^2)+(prob_change_in_y^2))/1000
    
    # Calculate the area based on a binary classification.
    test_bin <- BinaryTransformation(map, 0.5)
    prob_bin <- BinaryTransformation(prob_raster, 0.5)
    
    # Calculate range size change.
    prob_change <- as.data.frame(BIOMOD_RangeSize(prob_bin, test_bin)$Compt.By.Models)
    
    # Add the metrics to results.
    similarity_metrics <- data.frame(method = "BiasFile", 
                                     species = species,
                                     
                                     # Niche Similarity.
                                     p_schoeners = prob_overlaps$D,
                                     p_hellingers = prob_overlaps$I,
                                     p_spearmans = prob_overlaps$rank.cor,
                                     
                                     # Centroid shift.
                                     p_centroidshift = prob_magnitude,
                                     
                                     # Range size changes.
                                     p_loss = prob_change$PercLoss,
                                     p_gain = prob_change$PercGain,
                                     p_total_change = prob_change$SpeciesRangeChange,
                                     p_range_size = prob_change$CurrentRangeSize,
                                     
                                     stringsAsFactors=FALSE)
    
    write.csv(similarity_metrics, paste(file_pathway, "/", model_name, "_metrics.csv", sep=""), row.names = FALSE)
  }




                          ##### End of Script. #####
