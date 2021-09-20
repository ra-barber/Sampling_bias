################################################################################
                     ##### Cleaning occurrence data #####
################################################################################


# This script cleans occurrence data, removing duplicates for presence-only models.


# Load packages.
library(raster)
library(dplyr)
library(janitor)
library(stringr)




################################################################################
                     #### Load in hoverfly data ####


# Load all occurences from 1983 - 2002.
D1 <- read.csv("Data/Hoverfly Data/1983_1987.csv")
D2 <- read.csv("Data/Hoverfly Data/1988_1992.csv")
D3 <- read.csv("Data/Hoverfly Data/1993_1997.csv")
D4 <- read.csv("Data/Hoverfly Data/1998_2002.csv")

# Merge into one data frame.
hoverflies <- rbind(D1, D2, D3, D4)
head(hoverflies)
colnames(hoverflies) <- c("Species", "Easting", "Northing")

# Convert locations to meters matching OSGB.
hoverflies$Easting <- hoverflies$Easting * 1000
hoverflies$Northing <- hoverflies$Northing * 1000
head(hoverflies)

# Clean species names.
hoverflies$Species <- str_extract(hoverflies$Species, "^[A-Z][a-z]* [a-z]*") %>% str_replace(" ", "_")




################################################################################
     #### Remove duplicate occurrences for presence only modeling ####


# Check the number of duplicates.
duplicates <- get_dupes(hoverflies)

# Remove the duplicates.
hoverflies <- distinct(hoverflies)

# Double check.
duplicates <- get_dupes(hoverflies)




################################################################################
#### Remove any points that lie outside the range of environmental rasters ####


# Read in environmental rasters.
bioclim_filenames <- list.files("Data/Environmental Layers/", full.names = T)
bio_stack <- stack(lapply(bioclim_filenames, raster))

# Pull out all the coordinates of the climate stack that have occurrences.
hoverfly_points <- extract(bio_stack, hoverflies[,2:3])

# Get the row numbers of occurrences that have NA values for environmental layers.
row_numbers <- which(apply(is.na(hoverfly_points), 1,  sum) > 0)

# Remove the rows which don't fit environmental data.
cleaned_hoverflies <- hoverflies[-row_numbers,]

# Total all species and filter for those over 1000.
sum_hoverflies <- plyr::count(cleaned_hoverflies, "Species")
over_1000 <- filter(sum_hoverflies, freq > 1000)
print(over_1000)

# Export the occurence counts.
write.csv(over_1000, "Data/Hoverfly Data/hoverfly_occurrence_number.csv", row.names = F)

# Filter the data to just these species.
over_1000_data <- cleaned_hoverflies %>% filter(Species %in% over_1000$Species)

# Export the cleaned occurrences.
write.csv(cleaned_hoverflies, "Data/Hoverfly Data/D_CLEANED_1983_2002.csv", row.names = F)
write.csv(over_1000_data, "Data/Hoverfly Data/D_ALL_1983_2002.csv", row.names = F)




################################################################################
            #### Load in the sampling effort data to clean ####


# Loading sampling effort data.
sample_effort <- read.csv("Data/Sample Effort/year_1km_nVisits.csv")

# Filter to only years 83-02.
years83_02<- filter(sample_effort, year > 1982)
years83_02 <- filter(years83_02, year < 2003)

# Convert to OSGB.
years83_02$x <- years83_02$x * 1000
years83_02$y <- years83_02$y * 1000

# Pull out all the coordinates of the climate stack that have occurrences.
sample_points <- extract(bio_stack, years83_02[,3:4])

# Get the row numbers of occurrences that have NA values for environmental layers.
row_numbers <- which(apply(is.na(sample_points), 1,  sum) > 0)

# Remove the rows which don't fit environmental data.
cleaned_sample_effort<- years83_02[-row_numbers,]

# Export the cleaned occurrences.
write.csv(cleaned_sample_effort, "Data/Sample Effort/Sample_Effort_1983_2002.csv", row.names = F)




                      ##### End of Script. #####