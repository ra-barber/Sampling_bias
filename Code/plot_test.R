# Load packages for plotting maps
library(stringr)
library(rgdal)
library(raster)
library(rasterVis)
library(RColorBrewer)

# Set working directory
setwd("Sampling Bias")




#### Prepare data for plotting ####

# List all the pathways for the maps
map_filepathways <- list.files("Results/Model Maps", recursive = T)

# Read in all the rasters
for (x in map_filepathways){
  raster <- raster(paste("Results/Model Maps", x, sep="/"))
  name <- gsub(".tif", "", x)
  # Sort the rasters in logical order
  name <- str_replace(name, "_(?=true)", "_A_")
  name <- str_replace(name, "_(?=sample_effort_)", "_B_")
  name <- str_replace(name, "_(?=target)", "_C_")
  name <- str_replace(name, "_(?=Buffer)", "_D_")
  name <- str_replace(name, "_(?=pop)", "_E_")
  name <- str_replace(name, "_(?=travel)", "_F_")
  # Strip pathway information before assigning the name
  name <- str_remove(name, ".*/")
  assign(name, raster)
}

# Remove Pipizella_vidutata raster as it has too few training samples after models were ran
rm(Pipizella_viduata_A_true_bias_pred)
rm(Pipizella_viduata_B_sample_effort_pred)
rm(Pipizella_viduata_C_target_group_pred)
rm(Pipizella_viduata_D_Buffer_pred)
rm(Pipizella_viduata_E_pop_dens_pred)
rm(Pipizella_viduata_F_travel_pred)


# Create a list of rasters and a list of unique species
raster_names <- grep("pred", ls(), value = T)
species_names <- str_remove(raster_names, "_(?=[A-F]).*")
species_names <- unique(species_names)

# Create a vector of plot names for subtitles
plot_titles <- c("a) No Correction", "b) Sample Effort", "c) Target Group", "d) 10km Buffer", "e) Population Density", "f) Travel Time")

# Create a palette and then use color ramp so that it can be stretched by level plot
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

# Create the directory to store the maps
dir.create("Plots/Maps")



#### Plot the maps of each method in a single plot using levelplot ####

# Change lattice options to remove margins
lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)

# Loop through each species and create a plot of the six different bias correction methods
i = 1
for (x in species_names){
  # Grep the rasters from the species, and then stack them for level plot
  indices <- grep(x, raster_names)
  species_stack <- stack(eval(parse(text=raster_names[indices[1]])),
                         eval(parse(text=raster_names[indices[2]])),
                         eval(parse(text=raster_names[indices[3]])),
                         eval(parse(text=raster_names[indices[4]])),
                         eval(parse(text=raster_names[indices[5]])),
                         eval(parse(text=raster_names[indices[6]])))
  # Get the species name as the plot title
  title <- gsub("_", " ", x)
  # Use levelplot to create the plot
  myplot <- levelplot(species_stack, 
                      margin=F, 
                      colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.9), width=0.3),  # Sets the legend
                      axis.line=list(col='black'), # Gives the maps a black outline
                      scales=list(alternating=c(1,0,0)), # Get rid of coordinate labels
                      xlab = "Easting", ylab = "Northing",
                      col.regions=myPalette, # Add the colour
                      at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
                      main=list(label=title, cex=1, vjust=2), # Add the title
                      names.attr=plot_titles, # Add the plot subtitles for each map
                      par.strip.text=list(cex=0.9), # Make the subtitles smaller
                      layout=c(3,2))
  # Get rid of white space around the plot
  myplot$aspect.fill = T
  # Create the pathway name to save the file
  name <- paste("Figure_", i, "_SuppInfo_", sep="")
  i = i + 1
  filepathway <- paste("Plots/Maps/",name, x, ".tiff", sep = "")
  # Save the plots as tiffs
  tiff(filepathway, res=300, width=168, height=230, units="mm", compression = "lzw+p")
  print(myplot)
  dev.off()
} 

layout.scale.bar(myplot)
layout.north.arrow(type = 1)

myplot <- levelplot(species_stack, 
                    margin=F, 
                    colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.9), width=0.3),  # Sets the legend
                    axis.line=list(col='black'), # Gives the maps a black outline # Get rid of coordinate labels
                    col.regions=myPalette, # Add the colour
                    at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
                    main=list(label=title, cex=1, vjust=2), # Add the title
                    names.attr=plot_titles, # Add the plot subtitles for each map
                    par.strip.text=list(cex=0.9),
                    xscale.components = xscale.raster,
                    layout.scale.bar=list(height=0.05),# Make the subtitles smaller
                    layout=c(3,2))
levelplot(species_stack, 
          margin=F, 
          colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.9), width=0.3),  # Sets the legend
          axis.line=list(col='black'), # Gives the maps a black outline
          scales=list(draw=FALSE), # Get rid of coordinate labels
          col.regions=myPalette, # Add the colour
          at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
          main=list(label=title, cex=1, vjust=2), # Add the title
          names.attr=plot_titles, # Add the plot subtitles for each map
          par.strip.text=list(cex=0.9),
          auto.key=F, # Make the subtitles smaller
          layout=c(3,2))


spplot(species_stack,
       margin=F, 
       colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.9), width=0.3),  # Sets the legend
       axis.line=list(col='black'), # Gives the maps a black outline
       scales=list(draw=FALSE), # Get rid of coordinate labels
       col.regions=myPalette, # Add the colour
       at=seq(0,1, len=99), # Set the breaks (1 for each percentage)
       main=list(label=title, cex=1, vjust=2), # Add the title
       names.attr=plot_titles, # Add the plot subtitles for each map
       par.strip.text=list(cex=0.9),
       auto.key=F, # Make the subtitles smaller
       layout=c(3,2), 
       fill = NULL) + sp.theme(set=FALSE)

?sp.theme
  
  latticeExtra::layer(SpatialPolygonsRescale(layout.north.arrow(),
                                                                              offset = c(1200,20000), scale = 1000))
  
  
extent(species_stack)
dimensions(species_stack[[3]])

  latticeExtra::layer({
    SpatialPolygonsRescale(layout.north.arrow(),
                           offset = c(120000,20000))}, force=T, under = F)
  
  latticeExtra::layer({
            SpatialPolygonsRescale(layout.scale.bar(),
                                   offset = c(1200000,20000))}, packets = 1)
  
  levelplot(species_stack[[1]])  + latticeExtra::layer({
    SpatialPolygonsRescale(layout.north.arrow(),
                           offset = c(1200000,200000), scale=1000)}, style = T)
  
  devtools::install_github('oscarperpinan/rasterVis')
  
spatialra
?layer
?SpatialPolygonsRescale

trellis.par.set(layout.scale.bar())
?layer
?sp.panel.layout()
  layer({
                      xs <- seq(181000, 181400, by=100)
                      grid.rect(x=xs, y=330500,
                                width=100, height=30,
                                gp=gpar(fill=rep(c('transparent', 'black'), 2)),
                                default.units='native')
                      grid.text(x= xs - 50, y=330560, seq(0, 400, by=100),
                                gp=gpar(cex=1),
                                default.units='native')
                    })
myplot
?levelplot

myplot[1]
layout.scale.bar((species_stack[[6]]))

install.packages("ggsn")
library(ggsn)

levelplot(myplot, par.settings=layout.scale.bar())
library(sp)

#make a separate plot for use in the text
indices <- grep("Baccha", raster_names)
species_stack <- stack(eval(parse(text=raster_names[indices[1]])),
                       eval(parse(text=raster_names[indices[2]])),
                       eval(parse(text=raster_names[indices[3]])),
                       eval(parse(text=raster_names[indices[4]])),
                       eval(parse(text=raster_names[indices[5]])),
                       eval(parse(text=raster_names[indices[6]])))
title <- gsub("_", " ", x) #get the species name as the plot title
myplot <- levelplot(species_stack, 
                    margin=F, 
                    colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.6), width=0.3), #sets the legend
                    axis.line=list(col='black'), # give the maps a black outline
                    scales=list(draw=FALSE), #get rid of coordinate labels
                    col.regions=myPalette, #add the colour
                    at=seq(0,1, len=99), #set the breaks
                    names.attr=plot_titles, # add the plot titles for each map
                    par.strip.text=list(cex=0.6), # make the subtitles smaller
                    layout=c(3,2))
myplot$aspect.fill = T # get rid of white space around the plot
filepathway <- paste("Plots/Baccha_text_plot.tiff", sep = "")
tiff(filepathway, res=250, width=6, height=8, units="in")
print(myplot)
dev.off()




#plot rasters of the different bias files and sampling effort
sampling_effort <- raster("Bias Files/sample_effort_bias_file.tif")
population_density <- raster("Bias Files/pop_dens_bias_file.tif")
travel_time <- raster("Bias Files/travel_time_bias_file.tif")
target_group <- raster("Bias Files/Target Group/Baccha_elongata_target_group_bias_file.tif")
distance <- raster("Bias Files/Buffers/Baccha_elongata_buffer_biasfile.tif")

bias_files <- stack(sampling_effort, target_group, distance, population_density, travel_time)

bias_names <- c("a) Sampling Effort", "b) Target Group","c) 10km Buffer", "d) Population Density", "e) Travel Time")

myplot <- levelplot(bias_files, 
                    margin=F, 
                    colorkey=list(space="bottom", labels=list(at=0:1,font=3, cex=0.6), width=0.3), #sets the legend
                    axis.line=list(col='black'), # give the maps a black outline
                    scales=list(draw=FALSE), #get rid of coordinate labels
                    col.regions=myPalette, #add the colour
                    at=seq(0,1, len=99), #set the breaks
                    names.attr=bias_names, # add the plot titles for each map
                    par.strip.text=list(cex=0.6),
                    layout=c(3,2))

myplot$aspect.fill = T # get rid of white space around the plot
filepathway <- paste("Plots/biasfiles.tiff", sep = "")
tiff(filepathway, res=250, width=6, height=8, units="in")
print(myplot)
dev.off()