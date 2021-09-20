###############################################################################
                       ##### Functions #####

# 
#  spec <- readRDS("Data/Virtual/Species/random_species_19.rds")
# # 
#  spec
# # 
#  library(raster)

# spec2 <- readRDS("Data/Virtual_03_07_beta/Virtual/Species/random_species_24.rds")
# 
# spec2

# background_points <- randomPoints(mask = uniform_biasfile, n = 10000, prob = TRUE)
# 
# 
# res(buffer10k_stack[[1]])
# 
# sum(buffer10k_stack[[1]])

# Should check this to make sure it doesn't change when we set the seed and run 
# # not in parallel.
# cellStats(virtual_buffer10k_stack, "sum")

# For running maxent with virtual data.
virtual_maxent <- function(species, method, biasfile){
  
  
  ##### Modeling #####
  
  # Get the model name.
  model_name <- paste(species, method, sep="_")
  
  # Filter occurrences for species.
  #occurrences <- virtual_occurences %>% filter(random_name == species) %>% select(PO_x_biased, PO_y_biased)
  
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
  #map <- raster.transformation(map, trans = "norm")
  
  # Save the map.
  writeRaster(map, paste(file_pathway, "/", model_name, "_map.tif", sep=""), overwrite=TRUE)
  
  
  ##### Similarity Metrics ######
  
  # Get the reference raster.  #### Need to change to get the right species.
  ref_raster <- virtual_stack[[paste(species, "_Unbiased_map", sep="")]]
  #prob_raster <- probability_stack[[paste("prob_", species, sep="")]]
  
  # Calculate the niche overlap using Schoner's Distance (Package: ENMTools)
  overlaps <- raster.overlap(ref_raster, map)
  #prob_overlaps <- raster.overlap(prob_raster, map)
  
  # Calculate the centroids.
  ref_centroid <- wt.centroid(rasterToPoints(ref_raster, spatial = T))
  #prob_centroid <- wt.centroid(rasterToPoints(prob_raster, spatial = T))
  test_centroid <- wt.centroid(rasterToPoints(map, spatial = T))
  
  # Calculate differences.
  change_in_x <- xmin(test_centroid) - xmin(ref_centroid)
  change_in_y <- ymin(test_centroid) - ymin(ref_centroid)
  #prob_change_in_x <- xmin(test_centroid) - xmin(prob_centroid)
  #prob_change_in_y <- ymin(test_centroid) - ymin(prob_centroid)
  
  # Calculate the Euclidean distance in km.
  magnitude <- sqrt((change_in_x^2)+(change_in_y^2))/1000
  #prob_magnitude <- sqrt((prob_change_in_x^2)+(prob_change_in_y^2))/1000
  
  # Calculate the area based on a binary classification.
  ref_bin <- BinaryTransformation(ref_raster, 0.5)
  #prob_bin <- BinaryTransformation(prob_raster, 0.5)
  test_bin <- BinaryTransformation(map, 0.5)
  
  # Calculate range size change.
  ref_change <- as.data.frame(BIOMOD_RangeSize(ref_bin, test_bin)$Compt.By.Models)
  #prob_change <- as.data.frame(BIOMOD_RangeSize(prob_bin, test_bin)$Compt.By.Models)
  
  # # Find the threshold that minimizes gains and losses, so as to not bias either way.
  # thresholds <- seq(from=0.01, to=0.99, by=0.01)
  # 
  # # Create a data frame to store thresholds.
  # range_changes <- data.frame()
  # prob_range_changes <- data.frame()
  # 
  # # Loop through thresholds.
  # for (thresh in thresholds){
  #   # Calculate the area based on a binary classification.
  #   test_bin <- BinaryTransformation(map, thresh)
  #   
  #   # Calculate range size change.
  #   loop_change <- as.data.frame(BIOMOD_RangeSize(ref_bin, test_bin)$Compt.By.Models)
  #   prob_loop_change <- as.data.frame(BIOMOD_RangeSize(prob_bin, test_bin)$Compt.By.Models)
  #   
  #   loop_change <- cbind(thresh, loop_change)
  #   prob_loop_change <- cbind(thresh, prob_loop_change)
  #   
  #   range_changes <- rbind(range_changes, loop_change)
  #   prob_range_changes <- rbind(prob_range_changes, prob_loop_change)
  # }
  # 
  # # Pull out the row that minimises gains and losses. 
  # row_number <- which.min(abs(range_changes$PercLoss) + abs(range_changes$PercGain))
  # 
  # best_gain <- range_changes$PercGain[row_number]
  # best_loss <- range_changes$PercLoss[row_number]
  # best_change <- range_changes$SpeciesRangeChange[row_number]
  # best_threshold <- range_changes$thresh[row_number]
  # 
  # prob_row_number <- which.min(abs(prob_range_changes$PercLoss) + abs(prob_range_changes$PercGain))
  # 
  # prob_best_gain <- prob_range_changes$PercGain[prob_row_number]
  # prob_best_loss <- prob_range_changes$PercLoss[prob_row_number]
  # prob_best_change <- prob_range_changes$SpeciesRangeChange[prob_row_number]
  # prob_best_threshold <- prob_range_changes$thresh[prob_row_number]
  # 
  # 
  
  
  
  ###### Model Evaluation #####
  
  # # Get presences and absences.
  # virtual_occurences <- read.csv("Data/Virtual/Occurrences/all_occurrences.csv")
  # presences <- virtual_occurences %>% filter(random_name  == species & P_or_A_unbiased == 1) %>% select(PA_x_unbiased, PA_y_unbiased)
  # absences <- virtual_occurences %>% filter(random_name  == species & P_or_A_unbiased == 0) %>% select(PA_x_unbiased, PA_y_unbiased)
  # 
  # # Evaluate model using dismo.
  # model_eval <- dismo::evaluate(presences, absences, model, bio_stack)
  # 
  # # Get the threshold that maximises true postives and negatives.
  # threshold <- model_eval@t[which.max(model_eval@TPR + model_eval@TNR)]
  # confusion <- as.data.frame(cbind(model_eval@t, model_eval@confusion))
  # confusion <- confusion %>% filter(V1 == threshold)
  # 
  # # Get under and over-prediction metrics.
  # over_pred <- confusion$fp/(confusion$tp+confusion$fp)
  # under_pred <- confusion$fn/(confusion$tp+confusion$fn)
  # 
  # # Get sensitivity and specificity.
  # sensitivity <- confusion$tp/(confusion$tp+confusion$fn)
  # specifity <- confusion$tn/(confusion$tn+confusion$fp)
  # 
  # Add the metrics to results.
  similarity_metrics <- data.frame(method = method, 
                                   species = species,
                                   
                                   # Niche Similarity.
                                   schoeners = overlaps$D,
                                   hellingers = overlaps$I,
                                   spearmans = overlaps$rank.cor,
                                   # Niche Similarity.
                                   #p_schoeners = prob_overlaps$D,
                                   #p_hellingers = prob_overlaps$I,
                                   #p_spearmans = prob_overlaps$rank.cor,
                                   
                                   # Centroid shift.
                                   centroidshift = magnitude,
                                   #p_centroidshift = prob_magnitude,
                                   
                                   # Range size changes.
                                   loss = ref_change$PercLoss,
                                   gain = ref_change$PercGain,
                                   total_change = ref_change$SpeciesRangeChange,
                                   range_size = ref_change$CurrentRangeSize,
                                   
                                   # best_loss = best_loss,
                                   # best_gain = best_gain,
                                   # best_change = best_change,
                                   # best_threshold = best_threshold,
                                   # 
                                   #p_loss = prob_change$PercLoss,
                                   #p_gain = prob_change$PercGain,
                                   #p_total_change = prob_change$SpeciesRangeChange,
                                   #p_range_size = prob_change$CurrentRangeSize,
                                   # 
                                   # p_best_loss = prob_best_loss,
                                   # p_best_gain = prob_best_gain,
                                   # p_best_change = prob_best_change,
                                   # p_best_threshold = prob_best_threshold,
                                   # 
                                   # auc = model_eval@auc,
                                   # over_prediction = over_pred,
                                   # under_prediction = under_pred,
                                   # sensitivity = sensitivity,
                                   # specifity = specifity,
                                   # threshold = threshold,
                                   # 
                                   stringsAsFactors=FALSE)
  
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
  #map <- raster.transformation(map, trans = "norm")
  
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
  
  # # Find the threshold that minimizes gains and losses, so as to not bias either way.
  # thresholds <- seq(from=0.01, to=0.99, by=0.01)
  # 
  # # Create a data frame to store thresholds.
  # range_changes <- data.frame()
  # 
  # # Loop through thresholds.
  # for (thresh in thresholds){
  #   # Calculate the area based on a binary classification.
  #   test_bin <- BinaryTransformation(map, thresh)
  #   
  #   # Calculate range size change.
  #   loop_change <- as.data.frame(BIOMOD_RangeSize(ref_bin, test_bin)$Compt.By.Models)
  #   
  #   loop_change <- cbind(thresh, loop_change)
  #   range_changes <- rbind(range_changes, loop_change)
  # }
  # 
  # # Pull out the row that minimises gains and losses. 
  # row_number <- which.min(abs(range_changes$PercLoss) + abs(range_changes$PercGain))
  # 
  # best_gain <- range_changes$PercGain[row_number]
  # best_loss <- range_changes$PercLoss[row_number]
  # best_change <- range_changes$SpeciesRangeChange[row_number]
  # best_threshold <- range_changes$thresh[row_number]
  # 
  # Add the metrics to results.
  similarity_metrics <- data.frame(method = method, 
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
                                   
                                   # best_loss = best_loss,
                                   # best_gain = best_gain,
                                   # best_change = best_change,
                                   # best_threshold = best_threshold,
                                   
                                   stringsAsFactors=FALSE)
  
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