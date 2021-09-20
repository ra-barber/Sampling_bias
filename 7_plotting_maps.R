################################################################################
                        ##### Plotting maps #####
################################################################################


# This script generates all the maps using the predicted distributions from
# maxent models, as well as the bias file plots in the main text. The maps were
# also visually inspected using these plot settings.


# Load packages.
#library(stringr)
library(raster)
library(rasterVis)
library(RColorBrewer)




################################################################################
              ##### Read data and prepare for plots #####


##### Read in random species and supplementary info distributions ######

# Read in the probability rasters.
probability_filenames <- list.files("Data/Virtual/Probabilites/", full.names = T)
probability_stack <- stack(lapply(probability_filenames, raster))

# Read in suitability rasters.
suitability_filenames <- list.files("Data/Virtual/Suitability/", full.names = T)
suitability_stack <- stack(lapply(suitability_filenames, raster))

# Read in suitability rasters.
pa_filenames <- list.files("Data/Virtual/PA Rasters/", full.names = T)
pa_stack <- stack(lapply(pa_filenames, raster))

# Create a vector of folders that have virtual species maps used in supplementary.
reference_folders <- c("Results/Virtual/Biased", "Results/Virtual/BiasFile",
                       "Results/Virtual/Unbiased")

# Read in virtual species occurrences and create a vector of names to pull out species.
virtual_occurrences <- read.csv("Data/Virtual/Occurrences/all_occurrences.csv")
virtual_species <- unique(virtual_occurrences$random_name)
virtual_species <- paste(virtual_species, "_", sep="") # So that 1 doesn't pull out 10 etc.

# Read in reference results, using names of virtual species to grep from the whole list.
map_filenames <- list.files(reference_folders, recursive = T, full.names = T, pattern = ".tif")
reference_stacks <- lapply(virtual_species, function(x) grep(x, map_filenames, value = TRUE))
reference_stacks <- lapply(reference_stacks, stack)

# Add the probability and suitability stacks to the raster.
for (x in 1:50){
  name <- paste("suitab_random_species_", x, sep="")
  reference_stacks[[x]] <- stack(suitability_stack[[name]], reference_stacks[[x]])
  name <- paste("prob_random_species_", x, sep="")
  reference_stacks[[x]] <- stack(probability_stack[[name]], reference_stacks[[x]])
  name <- paste("pa_random_species_", x, sep="")
  reference_stacks[[x]] <- stack(pa_stack[[name]], reference_stacks[[x]])
}



##### Reading in virtual species ######

# Get the map files.
virtual_folders <- c("Results/Virtual/NoCorrection",
             "Results/Virtual/TargetGroup","Results/Virtual/Buffer10k",
             "Results/Virtual/PopDens", "Results/Virtual/TravelTime",
             "Results/Virtual/Unbiased")

map_filenames <- list.files(virtual_folders, recursive = T, full.names = T, pattern = ".tif")
virtual_stacks <- lapply(virtual_species, function(x) grep(x, map_filenames, value = TRUE))
virtual_stacks <- lapply(virtual_stacks, stack)

# # Add the probability and suitability stacks to the raster.
# for (x in 1:50){
#   name <- paste("prob_random_species_", x, sep="")
#   virtual_stacks[[x]] <- stack(probability_stack[[name]], virtual_stacks[[x]])
# }

##### Reading in empirical species ######

# Read in hoverflies.
hoverflies <- read.csv("Data/Hoverfly Data/D_ALL_1983_2002.csv")
hoverfly_species <- unique(hoverflies$Species)

# Get folders for maps.
folders <- c("Results/Emperical/SampleEffort", "Results/Emperical/NoCorrection",
             "Results/Emperical/TargetGroup","Results/Emperical/Buffer10k",
             "Results/Emperical/PopDens", "Results/Emperical/TravelTime")

# Read in the maps.
emperical_filenames <- list.files(folders, recursive = T, full.names = T, pattern = ".tif")
emperical_stacks <- lapply(hoverfly_species, function(x) grep(x, emperical_filenames, value = TRUE))
emperical_stacks <- lapply(emperical_stacks, stack)
names(virtual_stacks[[1]])


# Change the order of the stacks.
emperical_stacks <- lapply(emperical_stacks, function(x) stack(x[[4]], x[[2]], 
                                                               x[[5]], x[[1]],
                                                               x[[3]], x[[6]]))

virtual_stacks <- lapply(virtual_stacks, function(x) stack(x[[6]], x[[2]],
                                                           x[[4]], x[[1]],
                                                           x[[3]], x[[5]]))

reference_stacks <- lapply(reference_stacks, function(x) stack(x[[3]], x[[2]],
                                                               x[[1]], x[[4]],
                                                               x[[6]], x[[5]]))



################################################################################
                       ##### Plot all species #####


# Create a palette.
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

# Create plot subtitles.
plot_titles <- c("a) Sample Effort", "b) No Correction", "c) Target Group", 
                 "d) 10km Buffer", "e) Population Density", "f) Travel Time")

# Plot all empirical species on one PDF.
pdf("Plots/empirical_distributions.pdf")
for (x in 1:length(emperical_stacks)){
  title <- gsub("_", " ", hoverfly_species[x])
  plot <- levelplot(emperical_stacks[[x]], 
                    margin=F, 
                    colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.6), width=0.25),  # Sets the legend
                    axis.line=list(col='black'), # Gives the maps a black outline
                    scales=list(draw=FALSE), # Get rid of coordinate labels
                    col.regions=myPalette, # Add the colour
                    at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
                    main=list(label=title, cex=0.7, vjust=2), # Add the title
                    names.attr=plot_titles, # Add the plot subtitles for each map
                    par.strip.text=list(cex=0.6), # Make the subtitles smaller
                    layout=c(3,2)
  )
  print(plot, newpage = TRUE)
}
dev.off()


# library(tools) # Using compression via compact PDF doesn't seem to work.
# compactPDF("Plots/emperical_distributions.pdf", gs_quality = "screen", gs_cmd = "C:/Program Files/gs/gs9.54.0/bin/gswin64c.exe")

# Virtual Plot labels.
v_plot_titles <- c("a) Unbiased Occurrences", "b) No Correction", 
                   "c) Target Group", "d) 10km Buffer", 
                   "e) Population Density", "f) Travel Time")


# All as one pdf for virtual species.
pdf("Plots/virtual_distributions.pdf")
for (x in 1:length(virtual_stacks)){
  title <- paste("Random Species ", x, sep="")
  plot <- levelplot(virtual_stacks[[x]], 
                    margin=F, 
                    colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.6), width=0.25),  # Sets the legend
                    axis.line=list(col='black'), # Gives the maps a black outline
                    scales=list(draw=FALSE), # Get rid of coordinate labels
                    col.regions=myPalette, # Add the colour
                    at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
                    main=list(label=title, cex=0.7, vjust=2), # Add the title
                    names.attr=v_plot_titles, # Add the plot subtitles for each map
                    par.strip.text=list(cex=0.6), # Make the subtitles smaller
                    layout=c(3,2)
  )
  print(plot, newpage = TRUE)
}
dev.off()


# Random species plot labels.
r_plot_titles <- c("a) Suitability", "b) Probability", 
                   "c) Presence Absence", "d) Biased Dataset", 
                   "e) Unbiased Dataset", "f) Sample Effort")

# All as one pdf.
pdf("Plots/reference_distributions.pdf")
for (x in 1:length(reference_stacks)){
  title <- paste("Random Species ", x, sep="")
  plot <- levelplot(reference_stacks[[x]], 
                    margin=F, 
                    colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.6), width=0.25),  # Sets the legend
                    axis.line=list(col='black'), # Gives the maps a black outline
                    scales=list(draw=FALSE), # Get rid of coordinate labels
                    col.regions=myPalette, # Add the colour
                    at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
                    main=list(label=title, cex=0.7, vjust=2), # Add the title
                    names.attr=r_plot_titles, # Add the plot subtitles for each map
                    par.strip.text=list(cex=0.6), # Make the subtitles smaller
                    layout=c(3,2)
  )
  print(plot, newpage = TRUE)
}
dev.off()




################################################################################
              ##### Plots for species used in main text #####


# Change lattice options to remove margins.
lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)

# Plot empirical species for text.
myplot <- levelplot(emperical_stacks[[1]], 
                    margin=F, 
                    colorkey=list(space="right", labels=list(at=c(0.0, 0.5, 1.0),font=3, cex=1.2), width=0.9), #sets the legend
                    axis.line=list(col='black'), # give the maps a black outline
                    scales=list(draw=FALSE), #get rid of coordinate labels
                    col.regions=myPalette, #add the colour
                    at=seq(0,1, len=99), #set the breaks
                    names.attr=rep("", 6), # Remove plot labels for each map.
                    par.strip.text=list(cex=0.01), # make the subtitles smaller
                    par.settings = list(axis.line = list(col = "transparent"),
                                        strip.background = list(col = 'transparent'), 
                                        strip.border = list(col = 'transparent')),
                    layout=c(3,2))
myplot$aspect.fill = T # get rid of white space around the plot

# Save plot and add labels.
tiff("Plots/Baccha_text_plot.tiff", res=300, width=168, height=200, units="mm")
print(myplot)
trellis.focus("panel", 1, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(a)", cex = 1.5)
trellis.focus("panel", 2, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(b)", cex = 1.5)
trellis.focus("panel", 3, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(c)", cex = 1.5)
trellis.focus("panel", 1, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(d)", cex = 1.5)
trellis.focus("panel", 2, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(e)", cex = 1.5)
trellis.focus("panel", 3, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(f)", cex = 1.5)
trellis.unfocus()
dev.off()

# Plot virtual species for text.
myplot <- levelplot(virtual_stacks[[1]], 
                    margin=F, 
                    colorkey=list(space="right", labels=list(at=c(0.0, 0.5, 1.0),font=3, cex=1.2), width=0.9), #sets the legend
                    axis.line=list(col='black'), # give the maps a black outline
                    scales=list(draw=FALSE), #get rid of coordinate labels
                    col.regions=myPalette, #add the colour
                    at=seq(0,1, len=99), #set the breaks
                    names.attr=rep("", 6), # Remove plot labels for each map.
                    par.strip.text=list(cex=0.01), # make the subtitles smaller
                    par.settings = list(axis.line = list(col = "transparent"),
                                        strip.background = list(col = 'transparent'), 
                                        strip.border = list(col = 'transparent')),
                    layout=c(3,2))
myplot$aspect.fill = T # get rid of white space around the plot

# Save plot and add labels.
tiff("Plots/virtual_species_text_plot.tiff", res=300, width=168, height=200, units="mm")
print(myplot)
trellis.focus("panel", 1, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(a)", cex = 1.5)
trellis.focus("panel", 2, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(b)", cex = 1.5)
trellis.focus("panel", 3, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(c)", cex = 1.5)
trellis.focus("panel", 1, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(d)", cex = 1.5)
trellis.focus("panel", 2, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(e)", cex = 1.5)
trellis.focus("panel", 3, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(f)", cex = 1.5)
trellis.unfocus()
dev.off()




################################################################################
             ##### Plots of bias files used in main text #####


# Read in the bias files.
sampling_effort <- raster("Data/Bias Files/SampleEffort_biasfile.tif")
population_density <- raster("Data/Bias Files/PopDens_biasfile.tif")
travel_time <- raster("Data/Bias Files/TravelTime_biasfile.tif")
target_group <- raster("Data/Bias Files/TargetGroup_biasfile.tif")
distance <- raster("Data/Bias Files/Buffer 10k/Baccha_elongata_Buffer10k_biasfile.tif")
uniform <- raster("Data/Bias Files/Uniform_biasfile.tif")

# Stack the bias files in order.
bias_files <- stack(sampling_effort, uniform, target_group, distance, population_density, travel_time)

# Create a vector of labels.
# bias_names <- c("a) Uniform", 
#                 "b) Sampling Effort",
#                 "c) Target Group",
#                 "d) 10km Buffer", 
#                 "e) Population Density", 
#                 "f) Travel Time")

# Plot bias files.
myplot <- levelplot(bias_files, 
                    margin=F, 
                    colorkey=list(space="right", labels=list(at=c(0.0, 0.5, 1.0),font=3, cex=1.2), width=0.9), #sets the legend
                    axis.line=list(col='black'), # give the maps a black outline
                    scales=list(draw=FALSE), #get rid of coordinate labels
                    col.regions=myPalette, #add the colour
                    at=seq(0,1, len=99), #set the breaks
                    names.attr=rep("", 6), # Remove plot labels for each map.
                    par.strip.text=list(cex=0.01), # make the subtitles smaller
                    par.settings = list(axis.line = list(col = "transparent"),
                                        strip.background = list(col = 'transparent'), 
                                        strip.border = list(col = 'transparent')),
                    layout=c(3,2))
myplot$aspect.fill = T # get rid of white space around the plot

# Save plot and add labels.
tiff("Plots/biasfiles.tiff", res=300, width=168, height=200, units="mm")
print(myplot)
trellis.focus("panel", 1, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(a)", cex = 1.5)
trellis.focus("panel", 2, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(b)", cex = 1.5)
trellis.focus("panel", 3, 1, clip.off = TRUE)
panel.text(100000, 1200000, "(c)", cex = 1.5)
trellis.focus("panel", 1, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(d)", cex = 1.5)
trellis.focus("panel", 2, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(e)", cex = 1.5)
trellis.focus("panel", 3, 2, clip.off = TRUE)
panel.text(100000, 1200000, "(f)", cex = 1.5)
trellis.unfocus()
dev.off()




                       ##### End of Script. #####