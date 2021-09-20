################################################################################
                     ##### Statistical analysis #####
################################################################################


# This script performs all the statistical analysis and generates the figures 
# for those stats. Because our metrics were not normally distributed, we used
# a non-parametric friedman test, using repeated measures to account for the 
# effect of species and prevent pseudo-replication.


# Load packages.
library(rstatix)
library(ggpubr)
library(dplyr)




################################################################################
             ##### Read data and prepare for stats/plots #####


# Read in empirical similarity metrics.
results_filenames <- list.files("Results/Emperical", recursive = T, 
                                full.names = T, pattern = "_metrics.csv")
emperical_results <- bind_rows(lapply(results_filenames, read.csv))

# Create a list of folders used by virtual species in the main analysis.
virtual_folders <- c("Results/Virtual/NoCorrection/", "Results/Virtual/TargetGroup/",
                     "Results/Virtual/Buffer10k/", "Results/Virtual/PopDens/", 
                     "Results/Virtual/TravelTime/")

# Extract the virtual metrics data frames, and bind together.
results_filenames <- list.files(virtual_folders, recursive = T, 
                                full.names = T, pattern = "_metrics.csv")
virtual_results <- bind_rows(lapply(results_filenames, read.csv))

# Create list of folders used in supplementary results.
reference_folders <- c("Results/Virtual/Biased", "Results/Virtual/Unbiased/", 
                       "Results/Virtual/BiasFile/")

# Extract the results data frames and bind together.
reference_filenames <- list.files(reference_folders, recursive = T, 
                                  full.names = T, pattern = "_metrics.csv")
reference_results <- bind_rows(lapply(reference_filenames, read.csv))

# Convert strings to factors for stats.
emperical_results$species <- factor(emperical_results$species)
virtual_results$species <- factor(virtual_results$species)
reference_results$species <- factor(reference_results$species)

# Change the order of factors for plotting.
emperical_results$method <- factor(emperical_results$method, 
                                   levels = c("NoCorrection", "TargetGroup",
                                              "Buffer10k", "PopDens", 
                                              "TravelTime"))
virtual_results$method <- factor(virtual_results$method, 
                                 levels = c("NoCorrection", "TargetGroup",
                                            "Buffer10k", "PopDens", 
                                            "TravelTime"))
reference_results$method <- factor(reference_results$method, 
                                   levels = c("Biased", "Unbiased", "BiasFile"))

# Create plot labels.
plot_labels <- c("No\nCorrection", "Target\nGroup", 
                 "Radius", "Population\nDensity", "Travel\nTime")
ref_plot_labels <- c("Biased", "Unbiased", 
                     "Sample Effort")

# Save the plot theme.
plot_theme <- theme(axis.text.x=element_text(size=rel(0.85)), 
                    axis.text.y=element_text(size=rel(0.85)),
                    axis.title.y=element_text(size=rel(1)))




################################################################################
                    ##### Testing Empirical Species #####


# Do a friedman test of Schoener's distance.
e_schoeners_friedman <- emperical_results %>% friedman_test(schoeners ~ method | species)

# Record the effect size.
e_schoeners_effect_size <- emperical_results %>% friedman_effsize(schoeners ~ method | species)

# Generate pairwise comparisons.
e_schoeners_pair_wise <- emperical_results %>% wilcox_test(schoeners ~ method, paired = TRUE, p.adjust.method = "bonferroni")
e_schoeners_bias_wise <- e_schoeners_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
e_schoeners_plot <- ggboxplot(emperical_results, x = "method", y = "schoeners", add = "point", size = 0.4) +
  stat_pvalue_manual(e_schoeners_bias_wise, hide.ns = TRUE, y.position = c(0.98), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Schoener's Distance") + plot_theme


# Do a friedman test of centroid shift, with effect size and pairwise comparisons.
e_centroid_friedman <- emperical_results %>% friedman_test(centroidshift ~ method | species)
e_centroid_effect_size <- emperical_results %>% friedman_effsize(centroidshift ~ method | species)
e_centroid_pair_wise <- emperical_results %>% wilcox_test(centroidshift ~ method, paired = TRUE, p.adjust.method = "bonferroni")
e_centroid_bias_wise <- e_centroid_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
e_centroid_plot <- ggboxplot(emperical_results, x = "method", y = "centroidshift", add = "point", size = 0.4) +
  stat_pvalue_manual(e_centroid_bias_wise, hide.ns = TRUE, y.position = c(110), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Centroid Shift (km)") + plot_theme


# Do a friedman test of range gain, with effect size and pairwise comparisons.
e_gain_friedman <- emperical_results %>% friedman_test(gain ~ method | species)
e_gain_effect_size <- emperical_results %>% friedman_effsize(gain ~ method | species)
e_gain_pair_wise <- emperical_results %>% wilcox_test(gain ~ method, paired = TRUE, p.adjust.method = "bonferroni")
e_gain_bias_wise <- e_gain_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
e_gain_plot <- ggboxplot(emperical_results, x = "method", y = "gain", add = "point", size = 0.4) +
  stat_pvalue_manual(e_gain_bias_wise, hide.ns = TRUE, y.position = c(55), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Range Gain (%)") + plot_theme


# Do a friedman test of range loss, with effect size and pairwise comparisons.
e_loss_friedman <- emperical_results %>% friedman_test(loss ~ method | species)
e_loss_effect_size <- emperical_results %>% friedman_effsize(loss ~ method | species)
e_loss_pair_wise <- emperical_results %>% wilcox_test(loss ~ method, paired = TRUE, p.adjust.method = "bonferroni")
e_loss_bias_wise <- e_loss_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
e_loss_plot <- ggboxplot(emperical_results, x = "method", y = "loss", add = "point", size = 0.4) +
  stat_pvalue_manual(e_loss_bias_wise, hide.ns = TRUE, y.position = c(55), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Range Loss (%)") + plot_theme


# Plot the results together.
pdf("Plots/Statistics/emperical_results.pdf")
ggarrange(e_schoeners_plot + rremove("xlab"), e_centroid_plot + rremove("xlab"),
          e_gain_plot + rremove("xlab"), e_loss_plot + rremove("xlab"), 
          nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))
dev.off()

# Create seperate tiffs for publication.
tiff("Plots/emperical_results.tiff", res=900, width=200, height=250, units="mm")
ggarrange(e_schoeners_plot + rremove("xlab"), e_centroid_plot + rremove("xlab"),
          e_gain_plot + rremove("xlab"), e_loss_plot + rremove("xlab"), 
          nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))
dev.off()




################################################################################
                   ##### Testing Virtual Species #####


# Do a friedman test of Schoener's, with effect size and pairwise comparisons.
v_schoeners_friedman <- virtual_results %>% friedman_test(schoeners  ~ method | species)
v_schoeners_effect_size <- virtual_results %>% friedman_effsize(schoeners  ~ method | species)
v_schoeners_pair_wise <- virtual_results %>% wilcox_test(schoeners  ~ method, paired = TRUE, p.adjust.method = "bonferroni")
v_schoeners_bias_wise <- v_schoeners_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
v_schoeners_plot <- ggboxplot(virtual_results, x = "method", y = "schoeners", add = "point", size = 0.4) +
  stat_pvalue_manual(v_schoeners_bias_wise, hide.ns = TRUE, y.position = c(0.96), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Schoener's Distance") + plot_theme


# Do a friedman test of centroid, with effect size and pairwise comparisons.
v_centroidshift_friedman <- virtual_results %>% friedman_test(centroidshift ~ method | species)
v_centroidshift_effect_size <- virtual_results %>% friedman_effsize(centroidshift ~ method | species)
v_centroidshift_pair_wise <- virtual_results %>% wilcox_test(centroidshift ~ method, paired = TRUE, p.adjust.method = "bonferroni")
v_centroidshift_bias_wise <- v_centroidshift_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
v_centroidshift_plot <- ggboxplot(virtual_results, x = "method", y = "centroidshift", add = "point", size = 0.4) +
  stat_pvalue_manual(v_centroidshift_bias_wise, hide.ns = TRUE, y.position = c(150), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Centroid Shift (km)") + plot_theme


# Do a friedman test of range gain, with effect size and pairwise comparisons.
v_gain_friedman <- virtual_results %>% friedman_test(gain ~ method | species)
v_gain_effect_size <- virtual_results %>% friedman_effsize(gain ~ method | species)
v_gain_pair_wise <- virtual_results %>% wilcox_test(gain ~ method, paired = TRUE, p.adjust.method = "bonferroni")
v_gain_bias_wise <- v_gain_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
v_gain_plot <- ggboxplot(virtual_results, x = "method", y = "gain", add = "point", size = 0.4) +
  stat_pvalue_manual(v_gain_bias_wise, hide.ns = TRUE, y.position = c(60), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Range Gain (%)") + plot_theme  + coord_cartesian(ylim = c(0.05,60))


# Do a friedman test of range loss, with effect size and pairwise comparisons.
v_loss_friedman <- virtual_results %>% friedman_test(loss ~ method | species)
v_loss_effect_size <- virtual_results %>% friedman_effsize(loss ~ method | species)
v_loss_pair_wise <- virtual_results %>% wilcox_test(loss ~ method, paired = TRUE, p.adjust.method = "bonferroni")
v_loss_bias_wise <- v_loss_pair_wise %>% filter(group1 == "NoCorrection"|group2 == "NoCorrection")

# Create the plot.
v_loss_plot <- ggboxplot(virtual_results, x = "method", y = "loss", add = "point", size = 0.4) +
  stat_pvalue_manual(v_loss_bias_wise, hide.ns = TRUE, y.position = c(60), remove.bracket = TRUE, size = 5) +
  scale_x_discrete(labels=plot_labels) + 
  xlab("Correction Method") + ylab("Range Loss (%)") + plot_theme  + coord_cartesian(ylim = c(0.05,60))


# Plot results together.
pdf("Plots/Statistics/virtual_results.pdf")
ggarrange(v_schoeners_plot + rremove("xlab"), v_centroidshift_plot + rremove("xlab"),
          v_gain_plot + rremove("xlab"), v_loss_plot + rremove("xlab"), 
          nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))
dev.off()

tiff("Plots/virtual_results.tiff", res=900, width=200, height=250, units="mm")
ggarrange(v_schoeners_plot + rremove("xlab"), v_centroidshift_plot + rremove("xlab"),
          v_gain_plot + rremove("xlab"), v_loss_plot + rremove("xlab"), 
          nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))
dev.off()




################################################################################
                              # Export stats #

# Bind results together from the main analysis.
stats_results <- bind_rows(mget(ls(pattern = "friedman")))
effect_sizes <- bind_rows(mget(ls(pattern = "effect_size")))
statistic_results <- cbind(stats_results, effect_sizes)

# Export the results.
write.csv(statistic_results, "Results/main_statistics.csv", row.names = FALSE)

# Remove previous friedman and effect size objects before doing supplementary
# analysis.
rm(list=ls(pattern = "friedman"))
rm(list=ls(pattern = "effect_size"))




################################################################################
                    ##### Testing Sampling Effort #####


# Do a friedman test of Schoener's, with effect size and pairwise comparisons.
schoeners_friedman <- reference_results %>% friedman_test(p_schoeners ~ method | species)
schoeners_effect_size <- reference_results %>% friedman_effsize(p_schoeners ~ method | species)
schoeners_pair_wise <- reference_results %>% wilcox_test(p_schoeners ~ method, paired = TRUE, p.adjust.method = "bonferroni")

# Create the plot.
r_schoeners_plot <- ggboxplot(reference_results, x = "method", y = "p_schoeners", add = "point", size = 0.4) +
  stat_pvalue_manual(schoeners_pair_wise, hide.ns = TRUE, y.position = c(0.99, 1.02, 1.05), size = 5, bracket.size = 0.6) +
  scale_x_discrete(labels=ref_plot_labels) + 
  xlab("Correction Method") + ylab("Schoener's Distance") + plot_theme


# Do a friedman test of centroid shift, with effect size and pairwise comparisons.
r_centroid_friedman <- reference_results %>% friedman_test(p_centroidshift ~ method | species)
r_centroid_effect_size <- reference_results %>% friedman_effsize(p_centroidshift ~ method | species)
r_centroid_pair_wise <- reference_results %>% wilcox_test(p_centroidshift ~ method, paired = TRUE, p.adjust.method = "bonferroni")

# Create the plot.
r_centroid_plot <- ggboxplot(reference_results, x = "method", y = "p_centroidshift", add = "point", size = 0.4) +
  stat_pvalue_manual(r_centroid_pair_wise, hide.ns = TRUE, y.position = c(162, 174), size = 5, bracket.size = 0.6) +
  scale_x_discrete(labels=ref_plot_labels) + 
  xlab("Correction Method") + ylab("Centroid Shift") + plot_theme


# Do a friedman test of range gain, with effect size and pairwise comparisons.
r_gain_friedman <- reference_results %>% friedman_test(p_gain  ~ method | species)
r_gain_effect_size <- reference_results %>% friedman_effsize(p_gain  ~ method | species)
r_gain_pair_wise <- reference_results %>% wilcox_test(p_gain  ~ method, paired = TRUE, p.adjust.method = "bonferroni")

# Create the plot.
r_gain_plot <- ggboxplot(reference_results, x = "method", y = "p_gain", add = "point", size = 0.4) +
  stat_pvalue_manual(r_gain_pair_wise, hide.ns = TRUE, y.position = c(70), size = 5, bracket.size = 0.6) +
  scale_x_discrete(labels=ref_plot_labels) + 
  xlab("Correction Method") + ylab("Range Gain (%)") + plot_theme  + coord_cartesian(ylim = c(0.05,75))


# Do a friedman test of range loss, with effect size and pairwise comparisons.
r_loss_friedman <- reference_results %>% friedman_test(p_loss ~ method | species)
r_loss_effect_size <- reference_results %>% friedman_effsize(p_loss ~ method | species)
r_loss_pair_wise <- reference_results %>% wilcox_test(p_loss ~ method, paired = TRUE, p.adjust.method = "bonferroni")

# Create the plot.
r_loss_plot <- ggboxplot(reference_results, x = "method", y = "p_loss", add = "point", size = 0.4) +
  stat_pvalue_manual(r_loss_pair_wise, hide.ns = TRUE, y.position = c(63, 68, 73), size = 5, bracket.size = 0.6) +
  scale_x_discrete(labels=ref_plot_labels) + 
  xlab("Correction Method") + ylab("Range Loss (%)") + plot_theme  + coord_cartesian(ylim = c(0.05,75))


# Plot results together.
pdf("Plots/Statistics/reference_results.pdf", width = 5)
ggarrange(r_schoeners_plot + rremove("xlab"), r_centroid_plot + rremove("xlab"),
          r_gain_plot + rremove("xlab"), r_loss_plot + rremove("xlab"), 
          nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))
dev.off()

tiff("Plots/reference_results.tiff", res=900, width=200, height=250, units="mm")
ggarrange(r_schoeners_plot + rremove("xlab"), r_centroid_plot + rremove("xlab"),
          r_gain_plot + rremove("xlab"), r_loss_plot + rremove("xlab"), 
          nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))
dev.off()




################################################################################
                     # Export supplementary stats #

# Bind results together from the main analysis.
stats_results <- bind_rows(mget(ls(pattern = "friedman")))
effect_sizes <- bind_rows(mget(ls(pattern = "effect_size")))
statistic_results <- cbind(stats_results, effect_sizes)

# Export the results.
write.csv(statistic_results, "Results/supp_statistics.csv", row.names = FALSE)




                      ##### End of Script. #####