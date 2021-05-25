# Load packages needed for statistical analysis
library(tidyverse)
library(rstatix)
library(ggpubr)

# Set working directory
setwd("Sampling Bias")




#### Prepare data for stats/plotting ####

# Load in results
results <- read.csv("Results/Statistics/comparison_metrics.csv")

# Filter to remove sampling effort from analyses 
no_effort <- results %>% filter(method != "Sample Effort")


# Filter to remove models where the number of training samples dropped below 1000 due to NAs in environmental variables (these were coords outside model)
no_effort <- no_effort %>% filter(samplenumber > 999)

# Change centroid shift into km for easier interpretation
no_effort$centroidshift <- no_effort$centroidshift/1000
no_effort

# Plot histograms of variables
hist(no_effort$Schoeners)
hist(no_effort$rangesize)
hist(no_effort$centroidshift)

# Get a table of occurence number for each species
sample_number <- no_effort %>% filter(method == "True Bias")
sample_number <- sample_number[,c(2,4)]
sample_number$species <- gsub("_", " ", sample_number$species)
colnames(sample_number)[1] <- "Species"
colnames(sample_number)[2] <- "Occurences"
write.csv(sample_number, "Results/Statistics/occurence_number.csv", row.names = F)

# Convert strings to factors for stats
no_effort$method <- factor(no_effort$method)
no_effort$species <- factor(no_effort$species)

# Change the order of factors for plotting
no_effort$method <- factor(no_effort$method , levels=c("True Bias", "Target Group", "Buffer", "Population Density", "Travel Time"))

# Change scipen options for viewing signif values (optional)
options(scipen = 999)

# Check device is off
dev.off()

# Create a directory for storing plots 
dir.create("Plots")
dir.create("Plots/Statistics")




#### Do friedman tests of metrics and save results as a PDF ####

# Do a friedman test of schoeners
schoeners_friedman <- no_effort %>% friedman_test(Schoeners ~ method | species)
schoeners_friedman

# Look at the effect size
effect_size <- no_effort %>% friedman_effsize(Schoeners ~ method | species)
effect_size

# Look at the pairwise comparisons
pair_wise <- no_effort %>% wilcox_test(Schoeners ~ method, paired = TRUE, p.adjust.method = "bonferroni")
bias_wise <- pair_wise %>% filter(group1 == "True Bias"|group2 == "True Bias")

# Visulise as a box plot
pdf("Plots/Statistics/schoeners_boxplot.pdf")
bias_wise <- bias_wise %>% add_xy_position(x = "method")
ggboxplot(no_effort, x = "method", y = "Schoeners", add = "point", size = 0.5) +
  stat_pvalue_manual(bias_wise, hide.ns = TRUE, y.position = c(0.97, 0.985, 1.00, 1.015)) +
  scale_x_discrete(labels=c("No Correction", "Target Group",
                            "10km Buffer", "Population Density", "Travel Time")) + 
  xlab("Correction Method") + ylab("Schoener's Distance") + theme(axis.text.x=element_text(size=rel(0.75)))
dev.off()

# Now look at a change in range size
range_friedman <- no_effort %>% friedman_test(rangesize ~ method | species)
range_friedman

# Look at the effect size
effect_size <- no_effort %>% friedman_effsize(rangesize ~ method | species)
effect_size

# Look at the pairwise comparisons
pair_wise <- no_effort %>% wilcox_test(rangesize ~ method, paired = TRUE, p.adjust.method = "bonferroni")
bias_wise <- pair_wise %>% filter(group1 == "True Bias"|group2 == "True Bias")

# Visulise as a box plot
pdf("Plots/Statistics/rangesize_boxplot.pdf")
bias_wise <- bias_wise %>% add_xy_position(x = "method")
ggboxplot(no_effort, x = "method", y = "rangesize", add = "point", size = 0.5) +
  stat_pvalue_manual(bias_wise, hide.ns = TRUE, y.position = c(30, 35, 40)) +
  scale_x_discrete(labels=c("No Correction", "Target Group",
                            "10km Buffer", "Population Density", "Travel Time")) + 
  xlab("Correction Method") + ylab("Range Size Increase (%)") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + theme(axis.text.x=element_text(size=rel(0.75)))
dev.off()

# Look at a change in centroid shift
shift_friedman <- no_effort %>% friedman_test(centroidshift ~ method | species)
shift_friedman

# Look at the effect size
effect_size <- no_effort %>% friedman_effsize(centroidshift ~ method | species)
effect_size

# Look at the pairwise comparisons
pair_wise <- no_effort %>% wilcox_test(centroidshift ~ method, paired = TRUE, p.adjust.method = "bonferroni")
bias_wise <- pair_wise %>% filter(group1 == "True Bias"|group2 == "True Bias")

# Visulise as a box plot
pdf("Plots/Statistics/centroidshift_boxplot.pdf")
bias_wise <- bias_wise %>% add_xy_position(x = "method")
ggboxplot(no_effort, x = "method", y = "centroidshift", add = "point", size = 0.5) +
  stat_pvalue_manual(bias_wise, hide.ns = TRUE, y.position = c(125, 135, 145, 155), ) +
  scale_x_discrete(labels=c("No Correction", "Target Group",
                            "10km Buffer", "Population Density", "Travel Time")) + 
  xlab("Correction Method") + ylab("Centroid Shift (km)") + theme(axis.text.x=element_text(size=rel(0.75)))
dev.off()




#### Do one way Wilcoxon tests to see if results are significantly different from 1 or 0 ####

# Create object of methods to loop through
methods <- base::unique(no_effort$method)

# Create a dataframe for storing results of tests
wilcox_stats <- data.frame(dependent = character(),
                           method = character(),
                           statistic = numeric(),
                           p_value = numeric(),
                           ajusted_p <- numeric())

# Set levels of new dataframe factors
levels(wilcox_stats$dependent) <- c("Schoeners", "rangesize", "centroidshift")
levels(wilcox_stats$method) <- methods

# Compare to see if each method is significantly different from 1 (sample effort)
wilcox_stats[1:5,1] <- "Schoeners"
i=1
for (x in methods){
  # Do the one way wilcox text
  wilcox <- no_effort %>% filter(method==x) %>% wilcox_test(Schoeners ~ 1, mu = 1, alternative = "less", p.adjust.method = "bonferroni")
  # Save the method type
  wilcox_stats[i,2] <- x
  # Save the test statistic
  wilcox_stats[i,3] <- wilcox$statistic
  # Save the P value
  wilcox_stats[i,4] <- wilcox$p
  # Do a bonferroni correction to account for multiple tests
  new_p <- p.adjust(wilcox$p, method = "bonferroni", n=5)
  wilcox_stats[i,5] <- new_p
  i = i + 1
}

# Compare to see if each method is significantly different from 0 - (sample effort)
wilcox_stats[6:10,1] <- "rangesize"
for (x in methods){
  wilcox <- no_effort %>% filter(method==x) %>% wilcox_test(rangesize ~ 1, mu = 0, alternative = "two.sided", p.adjust.method = "bonferroni")
  wilcox_stats[i,2] <- x
  wilcox_stats[i,3] <- wilcox$statistic
  wilcox_stats[i,4] <- wilcox$p
  new_p <- p.adjust(wilcox$p, method = "bonferroni", n=5)
  wilcox_stats[i,5] <- new_p
  i = i + 1
}

# Compare to see if each method is significantly different from 0 - sample effort
wilcox_stats[11:15,1] <- "centroidshift"
for (x in methods){
  wilcox <- no_effort %>% filter(method==x) %>% wilcox_test(centroidshift ~ 1, mu = 0, alternative = "greater", p.adjust.method = "bonferroni")
  wilcox_stats[i,2] <- x
  wilcox_stats[i,3] <- wilcox$statistic
  wilcox_stats[i,4] <- wilcox$p
  new_p <- p.adjust(wilcox$p, method = "bonferroni")
  wilcox_stats[i,5] <- new_p
  i = i + 1
}

View(wilcox_stats)

