###############################################################################
                       ##### Functions #####

# For running maxent with virtual data.
virtual_maxent <- function(species, method, biasfile){
  
  
  ##### Modeling #####
  
  # Get the model name.
  model_name <- paste(species, method, sep="_")
  
  # Select 10000 background points without using probability.
  background_points <- randomPoints(mask = biasfile, n = 10000, prob = TRUE)
  
  # Create the file pathway.
  file_pathway <- paste("Results/Virtual", method, species, sep="/")
  
  # Create the maxent model with 10,000 background points
  model <- maxent(bio_stack, occurrences, background_points, factors="Land_Use", 
                  path=file_pathway, nbg=10000, overwrite=TRUE, args=c("maximumiterations=2000"))
  
  # Save the model as an RDS object.
  save(model, file = paste(file_pathway, "/", model_name, "_model.rds", sep=""))
  
  # Create the predicted map of each model.
  map <- predict(model, bio_stack, args=c('outputformat=cloglog'))
 
  # Save the map.
  writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  
  
  ##### Similarity Metrics ######
  
  # Get the reference raster.  #### Need to change to get the right species.
  ref_raster <- virtual_stack[[paste(species, "_Unbiased_map", sep="")]]
  
  # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
  overlaps <- raster.overlap(ref_raster, map)
  
  # Calculate the centroids.
  ref_centroid <- wt.centroid(rasterToPoints(ref_raster, spatial = T))
  test_centroid <- wt.centroid(rasterToPoints(map, spatial = T))
  
  # Calculate differences.
  change_in_x <- xmin(test_centroid) - xmin(ref_centroid)
  change_in_y <- ymin(test_centroid) - ymin(ref_centroid)
  
  # Calculate the Euclidean distance in km.
  magnitude <- sqrt((change_in_x^2)+(change_in_y^2))/1000
  
  # Calculate the area based on a binary classification.
  ref_bin <- BinaryTransformation(ref_raster, 0.5)
  test_bin <- BinaryTransformation(map, 0.5)
  
  # Calculate range size change.
  ref_change <- as.data.frame(BIOMOD_RangeSize(ref_bin, test_bin)$Compt.By.Models)

 
  # Add the metrics to results.
  similarity_metrics <- data.frame(# Iteration info.
                                   method = method, 
                                   species = species,
                                   
                                   # Niche Similarity.
                                   schoeners = overlaps$D,
                                   hellingers = overlaps$I,
                                   spearmans = overlaps$rank.cor,
                   
                                   # Centroid shift.
                                   centroidshift = magnitude,
                                   
                                   # Range size changes.
                                   loss = ref_change$PercLoss,
                                   gain = ref_change$PercGain,
                                   total_change = ref_change$SpeciesRangeChange,
                                   range_size = ref_change$CurrentRangeSize,
                                   
                                   stringsAsFactors=FALSE)
  
  # Export a csv.
  write.csv(similarity_metrics, paste(file_pathway, "/", model_name, "_metrics.csv", sep=""), row.names = FALSE)
}


# For running maxent with hoverfly data.
emperical_maxent <- function(species, method, biasfile){
  
  
  ##### Modeling #####
  
  # Get the model name.
  model_name <- paste(species, method, sep="_")
  
  # Filter occurrences for species.  #### Change selection.
  occurrences <- hoverflies %>% filter(Species == species) %>% select(Easting, Northing)
  
  # Select 10000 background points without using probability.
  background_points <- randomPoints(mask = biasfile, n = 10000, prob = TRUE)
  
  # Create the file pathway.
  file_pathway <- paste("Results/Emperical", method, species, sep="/")
  
  # Create the maxent model with 10,000 background points.
  model <- maxent(bio_stack, occurrences, background_points, factors="Land_Use", 
                  path=file_pathway, nbg=10000, overwrite=TRUE, args=c("maximumiterations=2000"))
  
  # Save the model as an RDS object.
  save(model, file = paste(file_pathway, "/", model_name, "_model.rds", sep=""))
  
  # Create the predicted map of each model
  map <- predict(model, bio_stack, args=c('outputformat=cloglog'))
  
  # Save the map.
  writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  
  ##### Similarity Metrics ######
  
  # Get the reference raster.  #### Need to change to get the right species.
  ref_raster <- emperical_stack[[paste(species, "_SampleEffort_map", sep="")]]
  
  # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
  overlaps <- raster.overlap(ref_raster, map)
  
  # Calculate the centroids.
  ref_centroid <- wt.centroid(rasterToPoints(ref_raster, spatial = T))
  test_centroid <- wt.centroid(rasterToPoints(map, spatial = T))
  
  # Calculate differences.
  change_in_x <- xmin(test_centroid) - xmin(ref_centroid)
  change_in_y <- ymin(test_centroid) - ymin(ref_centroid)
  
  # Calculate the Euclidean distance in km.
  magnitude <- sqrt((change_in_x^2)+(change_in_y^2))/1000
  
  # Calculate the area based on a binary classification.
  ref_bin <- BinaryTransformation(ref_raster, 0.5)
  test_bin <- BinaryTransformation(map, 0.5)
  
  # Calculate range size change.
  ref_change <- as.data.frame(BIOMOD_RangeSize(ref_bin, test_bin)$Compt.By.Models)
  
  # Add the metrics to results.
  similarity_metrics <- data.frame(# Iteration info.
                                   method = method, 
                                   species = species,
                                   
                                   # Niche Similarity.
                                   schoeners = overlaps$D,
                                   hellingers = overlaps$I,
                                   spearmans = overlaps$rank.cor,
                                   
                                   # Centroid shift.
                                   centroidshift = magnitude,
                                   
                                   # Range size changes.
                                   loss = ref_change$PercLoss,
                                   gain = ref_change$PercGain,
                                   total_change = ref_change$SpeciesRangeChange,
                                   range_size = ref_change$CurrentRangeSize,
                                   
                                   stringsAsFactors=FALSE)
  
  # Export a csv.
  write.csv(similarity_metrics, paste(file_pathway, "/", model_name, "_metrics.csv", sep=""), row.names = FALSE)
}



# For plotting colinearity. Edited from virtual species package.
remove_col_spearmans <- function (raster.stack, multicollinearity.cutoff = 0.7, select.variables = FALSE, 
                                  sample.points = FALSE, nb.points = 10000, plot = FALSE, 
                                  method = "pearson") 
{
  if (sample.points) {
    if (!is.numeric(nb.points)) {
      stop("nb.points must be a numeric value corresponding to the number of pixels to sample from raster.stack")
    }
    env.df <- sampleRandom(raster.stack, size = nb.points, 
                           na.rm = TRUE)
  }
  else {
    env.df <- getValues(raster.stack)
    if (any(is.na(env.df))) {
      env.df <- env.df[-unique(which(is.na(env.df), arr.ind = T)[, 
                                                                 1]), ]
    }
  }
  if (!is.numeric(multicollinearity.cutoff)) {
    stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")
  }
  else if (multicollinearity.cutoff > 1 | multicollinearity.cutoff < 
           0) {
    stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")
  }
  cor.matrix <- matrix(data = 0, nrow = nlayers(raster.stack), 
                       ncol = nlayers(raster.stack), dimnames = list(names(raster.stack), 
                                                                     names(raster.stack)))
  cor.matrix <- 1 - abs(stats::cor(env.df, method = method))
  dist.matrix <- stats::as.dist(cor.matrix)
  ahc <- stats::hclust(dist.matrix, method = "complete")
  groups <- stats::cutree(ahc, h = 1 - multicollinearity.cutoff)
  if (length(groups) == max(groups)) {
    message(paste("  - No multicollinearity detected in your data at threshold ", 
                  multicollinearity.cutoff, "\n", sep = ""))
    mc <- FALSE
  }
  else {
    mc <- TRUE
  }
  if (plot) {
    op <- par(no.readonly = TRUE)
    graphics::par(mar = c(5.1, 5.1, 4.1, 3.1))
    plot(ahc, hang = -1, xlab = "", ylab = "Distance (1 - Spearman's r)", 
         main = "", las = 1, sub = "", axes = F)
    graphics::axis(2, at = seq(0, 1, length = 6), las = 1)
    if (mc) {
      graphics::title(paste("Groups of intercorrelated variables at cutoff", 
                            multicollinearity.cutoff))
      par(xpd = T)
      rect.hclust(ahc, h = 1 - multicollinearity.cutoff)
    }
    else {
      graphics::title(paste("No intercorrelation among variables at cutoff", 
                            multicollinearity.cutoff))
    }
    par(op)
  }
  if (select.variables) {
    sel.vars <- NULL
    for (i in 1:max(groups)) {
      sel.vars <- c(sel.vars, sample(names(groups[groups == 
                                                    i]), 1))
    }
  }
  else {
    if (mc) {
      sel.vars <- list()
      for (i in groups) {
        sel.vars[[i]] <- names(groups)[groups == i]
      }
    }
    else {
      sel.vars <- names(raster.stack)
    }
  }
  return(sel.vars)
}
