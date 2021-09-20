################################################################################
                    ##### Creating virtual species #####
################################################################################


# This script creates virtual species for use in maxent models. All the virtual
# species were created using the environmental data produced in the 1st script. 
# Species were generated randomly, with a species prevalence of 0.5. This script
# uses parallel computing to improve processing times. It is possible to run 
# for each loops sequentially, by not setting up the number of cores using 
# the registerDoParallel() function. This may be necessary for some machines.


# Load packages.
library(raster)
library(dplyr)
library(ks)
library(spatialEco)
library(virtualspecies)
library(foreach)
library(doParallel)
library(biomod2)




################################################################################
                        ##### Read in Data #######


# Read in environmental rasters.
bioclim_filenames <- list.files("Data/Environmental Layers/", full.names = T)
bio_stack <- stack(lapply(bioclim_filenames, raster))

# Load sampling effort data.
sample_effort <- read.csv("Data/Sample Effort/Sample_Effort_1983_2002.csv")

# Set the seed. This doesn't create identical virtual species if they are
# generated in parallel.
set.seed(12345)




################################################################################
              ##### Creating sampling effort bias file #######


# Aggregate visits across different years for each grid cell.
sum_effort <- aggregate(sample_effort$nVisits, sample_effort[,3:4], FUN = sum)

# Create a matrix of coordinates.
coords <- cbind(sum_effort[,1], sum_effort[,2])

# Scale weights so that sum of all sampling weights = total sample size.
scale <- length(sum_effort[,3]) / sum(sum_effort[,3])
sum_effort[,4] <- sum_effort[,3] * scale

# Change names of columns.
colnames(sum_effort) <- c("Northing", "Easting", "nVisits", "Weights")
head(sum_effort)

# Do a 2d kernel density estimate (using scaled sampling effort as weighting) 
# to create a probability surface of sampling effort.
sample_density <- kde(coords, w=sum_effort$Weights)

# Create a raster of the sampling effort surface.
sample_raster <- raster(sample_density)
plot(sample_raster)

# Define sampling effort raster as in OSGB.
crs(sample_raster) <- '+init=EPSG:27700'

# Clip sample effort to same extent and resolution as land use data.
sample_raster <- resample(sample_raster, bio_stack, method='bilinear')

# Mask sample effort to remove outside area.
for (x in 1:nlayers(bio_stack)){
  sample_raster <- mask(sample_raster, bio_stack[[x]])
}

# Normalize bias file.
sample_raster <- raster.transformation(sample_raster, trans="norm")

# Check that the plots match up.
plot(bio_stack[[1]])
plot(sample_raster)
sample_raster

# Create a folder to store the bias files.
if (!dir.exists("Data/Bias Files")){
  dir.create("Data/Bias Files", recursive = TRUE)
}

# Export the bias file.
writeRaster(sample_raster, "Data/Bias Files/SampleEffort_biasfile.tif", overwrite=TRUE)

# You can skip the first step if the bias file already exists. 
# sample_raster <- raster("Data/Bias Files/SampleEffort_biasfile.tif")




# ################################################################################
             ##### Generate random species & occurrences #######


# Create a folder to store the random species, occurrences and rasters.
if (!dir.exists("Data/Virtual")){
  dir.create("Data/Virtual/Species", recursive = TRUE)
  dir.create("Data/Virtual/Occurrences", recursive = TRUE)
  dir.create("Data/Virtual/Suitability", recursive = TRUE)
  dir.create("Data/Virtual/Probabilites", recursive = TRUE)
  dir.create("Data/Virtual/PA Rasters", recursive = TRUE)
}

# Set up cluster for parallel processing.
#registerDoParallel(cores = 8)

# Data frames for storing occurrences and random species parameters.
virtual_occurences <- data.frame()
random_species_parameters <- data.frame()

# Loop through and create 50 random species. Save rasters and extract occurrences.
# virtual_occurences <- foreach(x=1:50,
#                               .inorder = FALSE,
#                               .errorhandling = "stop",
#                               .packages = c("raster", "virtualspecies", "biomod2"),
#                               .combine = rbind) %dopar% 
#   {

for (x in 1:50){
    # Create a name for the random species.
    random_name <- paste("random_species_", x, sep="")
    
    # Make sure the random species predicted suitability is not excessively small 
    # as this leads to models that perform poorly regardless of sampling bias, 
    # and is a poor validation of hoverfly species. Rerun random species if less
    # than 30% of cells (75000) have a suitability value below 0.5.
    # area <- 0
    # while(area < 75000){
    #   # Generate a random species object.
    #   random_species <- generateRandomSp(bio_stack, 
    #                                    species.prevalence = 0.5, 
    #                                    plot = FALSE, 
    #                                    niche.breadth = "wide")
    #   binary_map <- BinaryTransformation(random_species$suitab.raster, 0.5)
    #   area <- sum(na.omit(values(binary_map)))
    #  }
    
    # Try it this way.
    beta_value <- runif(1, min=0.3, max=0.7)
    random_species <- generateRandomSp(bio_stack,
                                       beta = beta_value,
                                       plot = FALSE, 
                                       niche.breadth = "wide")
    
    # Save the object.
    saveRDS(random_species, paste("Data/Virtual/Species/", random_name, ".rds", sep=""))
    
    # Save the parameter values.
    PA_beta <- random_species$PA.conversion[[4]]
    PA_prevalence <- random_species$PA.conversion[[5]]
    
    # Save the random species parameters used.
    loop_parameters <- cbind(random_name, PA_beta, PA_prevalence)
    random_species_parameters <- rbind(random_species_parameters, loop_parameters)
    
    # Save the random species rasters.
    writeRaster(random_species$suitab.raster, paste("Data/Virtual/Suitability/suitab_", random_name, ".tif", sep=""), overwrite=TRUE)
    writeRaster(random_species$probability.of.occurrence, paste("Data/Virtual/Probabilites/prob_", random_name, ".tif", sep=""), overwrite=TRUE)
    writeRaster(random_species$pa.raster, paste("Data/Virtual/PA Rasters/pa_", random_name, ".tif", sep=""), overwrite=TRUE)
    
    # Sample unbiased presence only points.
    PO_points <- sampleOccurrences(random_species, n = 1000, type = "presence only", plot = FALSE)
    
    # Sample unbiased presence-absence data we can use to test bias corrections.
    # sample.prevalence means we sample 500 presences and 500 absences.
    # PA_points <- sampleOccurrences(random_species, n = 1000, type = "presence-absence", 
    #                                sample.prevalence = 0.5, plot = FALSE)
    
    # Sample presence only points biased by sampling effort to create a biased dataset to test.
    PO_bias_points <- sampleOccurrences(random_species, n = 1000, type = "presence only",
                                        bias = "manual", weights = sample_raster, plot = FALSE)
    
    # Retrieve coordinates and change names.
    PO_points <- PO_points$sample.points[,1:2]
    colnames(PO_points) <- c("PO_x_unbiased", "PO_y_unbiased")
    
    PO_bias_points <- PO_bias_points$sample.points[,1:2]
    colnames(PO_bias_points) <- c("PO_x_biased", "PO_y_biased")
    
    # Should remove because we don't use them.
    #PA_points <- PA_points$sample.points[,1:3]
    #colnames(PA_points) <- c("PA_x_unbiased", "PA_y_unbiased", "P_or_A_unbiased")
    
    # Combine occurrences and save each species while looping.
    single_occurences <- cbind(random_name, PO_points, PO_bias_points)
    file_pathway <- paste("Data/Virtual/Occurrences/occurrences_", 
                          random_name, ".csv", sep="")
    write.csv(single_occurences, file_pathway, row.names = FALSE)
    
    # Return the data frame.
    virtual_occurences <- rbind(virtual_occurences, single_occurences)
    random_species_parameters <- rbind(random_species_parameters)
    
  }

# Save all the virtual species.
write.csv(virtual_occurences, "Data/Virtual/Occurrences/all_occurrences.csv", row.names = FALSE)

write.csv(random_species_parameters, "Data/Virtual/Occurrences/random_species_parameters.csv", row.names = FALSE)



                        ##### End of Script. #####

# for (x in 1:50){
#   binary_map <- BinaryTransformation(suitability_stack[[x]], 0.5)
#   area <- sum(na.omit(values(binary_map)))
#   print(area)
# }
# binary_map <- BinaryTransformation(random_species$suitab.raster, 0.5)
# area <- sum(na.omit(values(binary_map)))
# 
# 
# suitability_filenames <- list.files("Data/Virtual/Suitability/", full.names = T)
# suitability_stack <- stack(lapply(suitability_filenames, raster))
# plot(suitability_stack)
# plot(random_species_parameters$PA_beta, random_species_parameters$PA_prevalence)
# 
