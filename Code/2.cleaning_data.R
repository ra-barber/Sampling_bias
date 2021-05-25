# Load packages needed for spatial analysis
library(rgdal)
library(raster)
library(rgeos)


# Set working directory to the sampling bias folder
setwd("Sampling Bias")




#### Load in Hoverfly Data ####

# Load all occurences from 1983 - 2002
D3 <- read.csv("Data/Hoverfly Data/D3_1983_1987.csv")
D4 <- read.csv("Data/Hoverfly Data/D4_1988_1992.csv")
D5 <- read.csv("Data/Hoverfly Data/D5_1993_1997.csv")
D6 <- read.csv("Data/Hoverfly Data/D6_1998_2002.csv")


# Merge into one data.frame
hoverflies <- rbind(D3, D4, D5, D6)
head(hoverflies)

# Convert locations to metres
hoverflies$xx <- hoverflies$xx * 1000
hoverflies$yy <- hoverflies$yy * 1000
head(hoverflies)

# Remove duplicate ocurrences for presence only modelling
library(dplyr)
hoverflies <- distinct(hoverflies)




#### Remove any points that lie outside the range of environmental rasters ####
land_use <- raster("Data/Land Use/landuse.tif")

# Load climate file names
climate_filenames <- paste("Data/World_clim/OSGB/", list.files("Data/World_clim/OSGB/"), sep="")

# Read in climate rasters
for (x in 1:length(climate_filenames)){
  name <- paste("world_clim", x, sep="")
  raster <- raster(climate_filenames[x])
  assign(name, raster)
}

# Stack cliamte rasters
climate_stack <- stack(land_use, world_clim1, world_clim2, world_clim3, world_clim4,
                       world_clim5, world_clim6, world_clim7, world_clim8, 
                       world_clim9, world_clim10, world_clim11, world_clim12, 
                       world_clim13, world_clim14, world_clim15, world_clim16, 
                       world_clim17, world_clim18, world_clim19)



# Pull out all the coordinates of the climate stack that have occurrences
hoverfly_points <- extract(climate_stack, hoverflies[,2:3])

# Get the rownumbers of occurrences that have NA values for environmental layers
row_numbers <- which(apply(is.na(hoverfly_points), 1,  sum) > 0)

# Remove the rows which don't fit environmental data
cleaned_hoverflies <- hoverflies[-row_numbers,]

# Export the cleaned occurrences 
write.csv(cleaned_hoverflies, "Data/Hoverfly Data/D_ALL_1983_2002.csv", row.names = F)




#### Load in the sampling effort data to clean ####

# Loading sampling effort data
Sample_effort <- read.csv("Data/Sample_Effort/year_1km_nVisits.csv")

# Filter to only years 83-02
years83_02<- filter(Sample_effort, year > 1982)
years83_02 <- filter(years83_02, year < 2003)

# Convert to OSGB 
years83_02$x <- years83_02$x * 1000
years83_02$y <- years83_02$y * 1000

# Pull out all the coordinates of the climate stack that have occurrences
sample_points <- extract(climate_stack, years83_02[,3:4])

# Get the rownumbers of occurrences that have NA values for environmental layers
row_numbers <- which(apply(is.na(sample_points), 1,  sum) > 0)

# Remove the rows which don't fit environmental data
cleaned_sample_effort<- years83_02[-row_numbers,]

# Export the cleaned occurrences 
write.csv(cleaned_sample_effort, "Data/Sample_Effort/Sample_Effort_1983_2002.csv", row.names = F)



