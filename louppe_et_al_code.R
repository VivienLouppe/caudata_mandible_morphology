##### data preparation ============================================================================

### phylogeny

library(ape)
phy <- read.nexus(file = "./phylogeny_data.nexus")



### data

load(file = "./data_all.Rdata")
load(file = "./data_dent.Rdata")
load(file = "./data_artang.Rdata")







##### fitting evolutionary models ==================================================================

# Installing mvMORPH
library(devtools)
install_github("JClavel/mvMORPH", build_vignettes = FALSE)
library(mvMORPH)



### identify the best model

# Brownian motion
fit_bm <- mvgls(Y ~ lc * csize,
                data = data_all,
                tree = phy,
                model = "BM",
                method = "PL-LOOCV")
# Ornstein-Uhlenbeck 
fit_ou <- mvgls(Y ~ lc * csize,
                data = data_all,
                tree = phy,
                model = "OU",
                method = "PL-LOOCV")
# early burst 
fit_eb <- mvgls(Y ~ lc * csize,
                data = data_all,
                tree = phy,
                model = "EB",
                method = "PL-LOOCV")

# Pagel's lambda - a measure of “phylogenetic signal
fit_pl <- mvgls(Y ~ lc * csize,
                data = data_all,
                tree = phy,
                model = "lambda",
                method = "PL-LOOCV")

# Compare the models using the GIC criterion 
GIC(fit_bm)
GIC(fit_ou)
GIC(fit_eb)
GIC(fit_pl)



### fit complete models

# fit
fit_lc_size_all <- mvgls(Y ~ lc * csize,
                         data = data_all,
                         tree = phy,
                         model = "lambda",
                         method = "PL-LOOCV")
fit_lc_size_dent <- mvgls(Y ~ lc * csize,
                          data = data_dent,
                          tree = phy,
                          model = "lambda",
                          method = "PL-LOOCV")
fit_lc_size_artang <- mvgls(Y ~ lc * csize,
                            data = data_artang,
                            tree = phy,
                            model = "lambda",
                            method = "PL-LOOCV")


# list
fit_list_grp_size <- list(fit_lc_size_all = fit_lc_size_all,
                          fit_lc_size_dent = fit_lc_size_dent,
                          fit_lc_size_artang = fit_lc_size_artang)





##### Phylogenetic MANOVA and MANCOVA ==========================================================================

### mancovas

# create a list to store the results of the analyses
fit_list <- fit_list_grp_size
names_fit_list <- names(fit_list)
mancova_results <- list()

# loop to compute the analyses
for (i in 1:length(fit_list)) {
  model <- fit_list[[i]]
  result <- manova.gls(model,
                       test = "Wilks",
                       nperm = 999,
                       verbose = TRUE,
                       nbcores = 4L, # to speed up the calculation
                       type = "II")
  
  # store the results
  mancova_results[[names_fit_list[i]]] <- result
}



### test for effect size

# create a list to store the results of the analyses
names_fit_list <- names(fit_list)
effectsize_mancova <- list()

# loop to compute the analyses
for (i in 1:length(mancova_results)) {
  model <- mancova_results[[i]]
  result <- effectsize(model)
  
  # store the results
  effectsize_mancova[[names_fit_list[i]]] <- result
}







##### test for allometry =======================================================================================

library(mvMORPH)

### prepare data

# create a list to store the results of the analyses
load(file = "./data/mandibles/models_fits/fit_list_grp_size.Rdata")
fit_list <- list(fit_lc_size_all = fit_list_grp_size$fit_lc_size_all,
                 fit_lc_size_dent = fit_list_grp_size$fit_lc_size_dent,
                 fit_lc_size_artang = fit_list_grp_size$fit_lc_size_artang)
names_fit_list <- names(fit_list)



### contrast coding

# IMPORTANT: check the order of the coefficient to be sure to look at the interaction between the group and size
# allows to test hypothesis of allometry in the different groups
allom_bi <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0), nrow = 1) 
allom_dd <- matrix(c(0, 0, 0, 0, 0, 1, 0, 0), nrow = 1)
allom_pd <- matrix(c(0, 0, 0, 0, 0, 0, 1, 0), nrow = 1) 
allom_vi <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1), nrow = 1) 

contrast_list <- list(
  allom_bi = allom_bi,
  allom_dd = allom_dd,
  allom_pd = allom_pd,
  allom_vi = allom_vi
)
names_contrast <- names(contrast_list)



### function to test allometry

test_allometry <- function(model, contrast) {
  manova.gls(model,
             test = "Wilks",
             nperm = 999,
             L = contrast,
             verbose = TRUE,
             nbcores = 4L, # to speed up the calculation
             type = "II")
}



### run analysis

allometry_results <- list()

# loop to compute the analyses
for (i in 1:length(fit_list)) {
  for (j in 1:length(contrast_list)) {
    result <- test_allometry(fit_list[[i]], contrast_list[[j]])
    
    # store the results
    allometry_results[[paste(names_fit_list[i], names_contrast[j], sep = "_")]] <- result
  }
  
}



### test for effect size

load(file = "./results/allometry/allometry_results.Rdata")

# create a list to store the results of the analyses
names_allometry_results <- names(allometry_results)
effectsize_allom <- list()

# loop to compute the analyses
for (i in 1:length(allometry_results)) {
  model <- allometry_results[[i]]
  result <- effectsize(model)
  
  # store the results
  effectsize_allom[[names_allometry_results[i]]] <- result
}

capture.output(effectsize_allom, 
               file = "./results/allometry/effect_size_allometry.txt")



### extract the results 

allometry_results_table <- data.frame(
  Model = names(allometry_results),                             
  Test_stat = sapply(allometry_results, function(x) x$stat), 
  Effect_size = sapply(effectsize_allom, function(x) x$effect),
  Pvalue = sapply(allometry_results, function(x) x$pvalue),
  row.names = NULL
)






##### dfa =====================================================================================================

### compute dfa

# loop to compute the analyses
fit_list <- fit_list_grp_size
names_fit_list <- names(fit_list)
dfa_results <- list()

library(mvMORPH)
for (i in 1:length(fit_list)) {
  model <- fit_list[[i]]
  result <- mvgls.dfa(model)
  
  # store the results
  dfa_results[[names_fit_list[i]]] <- result
}



### extract results

dfa_scores <- list()
data_dfa <- list()
caudata_data <- read.csv2(file = "./data/mandibles/ident_tables/ident_ordered.csv")
d <- caudata_data[, -c(1, 3, 4, 7)]

for (i in 1:length(fit_list)) {
  # extract dfa scores
  dfa_scores[[names_fit_list[i]]] <- as.data.frame(dfa_results[[i]]$scores)
  
  # get rownames
  sp <- rownames(dfa_scores[[names_fit_list[i]]])
  
  # assemble scores and data
  data_dfa[[names_fit_list[i]]] <- cbind(d, dfa_scores[[i]])
  
  # remove rownames
  rownames(data_dfa[[names_fit_list[i]]]) <- NULL
}



### plot

data_dfa_lc <- list(all = data_dfa$fit_lc_size_all,
                    dent = data_dfa$fit_lc_size_dent,
                    artang = data_dfa$fit_lc_size_artang)
dfa_plots <- list()

for (i in seq_along(data_dfa_lc)) {
  # set theme
  theme_set(theme_bw())
  
  d <- data_dfa_lc[[i]]
  
  # plot df1-2
  plot1 <- ggplot(data = d,
                  aes(x = `Discr. 1`,
                      y = `Discr. 2`,
                      color = lc,
                      shape = family)) +
    geom_point(size = 1.1,
               stroke = 0.25,
               alpha = 0.7) + 
    # add data on axes
    geom_rug(show.legend = FALSE) + 
    
    # change legend color and text
    scale_color_manual(name = "Life cycle",
                       breaks = c("bi", "dd", "pd", "vi"),
                       labels = c("biphasic", "direct development", "paedomorphic", "viviparous"),
                       values = c("#d41159", "#ffc20a", "#1a85ff", "black")) +
    scale_shape_manual(name = "Family",
                       breaks = c("Ambystomatidae",
                                  "Amphiumidae", 
                                  "Cryptobranchidae", 
                                  "Dicamptodontidae", 
                                  "Hynobiidae", 
                                  "Plethodontidae", 
                                  "Proteidae", 
                                  "Rhyacotritonidae", 
                                  "Salamandridae", 
                                  "Sirenidae"),
                       labels = c("Ambystomatidae",
                                  "Amphiumidae", 
                                  "Cryptobranchidae", 
                                  "Dicamptodontidae", 
                                  "Hynobiidae", 
                                  "Plethodontidae", 
                                  "Proteidae", 
                                  "Rhyacotritonidae", 
                                  "Salamandridae", 
                                  "Sirenidae"),
                       values = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)) +
    
    # add title
    labs(
      x = "Discriminant 1",
      y = "Discriminant 2") + 
    
    # custom theme
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.caption = element_text(hjust = 0),
          legend.key = element_rect(fill = NA),
          legend.position = "none", # remove legend
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          # plot.margin = margin(
          #   t = 0.5,
          #   r = 0.5,
          #   b = 0.5,
          #   l = 0.5,
          #   unit = "cm"
          # ),
          aspect.ratio = 1)
  
  # plot df1-3
  plot2 <- ggplot(data = d,
                  aes(x = `Discr. 1`,
                      y = `Discr. 3`,
                      color = lc,
                      shape = family)) +
    geom_point(size = 1.1,
               stroke = 0.25,
               alpha = 0.7) +
    # add data on axes
    geom_rug(show.legend = FALSE) +
    
    # change legend color and text
    scale_color_manual(name = "Life cycle",
                       breaks = c("bi", "dd", "pd", "vi"),
                       labels = c("biphasic", "direct development", "paedomorphic", "viviparous"),
                       values = c("#d41159", "#ffc20a", "#1a85ff", "black")) +
    scale_shape_manual(name = "Family",
                       breaks = c("Ambystomatidae",
                                  "Amphiumidae",
                                  "Cryptobranchidae",
                                  "Dicamptodontidae",
                                  "Hynobiidae",
                                  "Plethodontidae",
                                  "Proteidae",
                                  "Rhyacotritonidae",
                                  "Salamandridae",
                                  "Sirenidae"),
                       labels = c("Ambystomatidae",
                                  "Amphiumidae",
                                  "Cryptobranchidae",
                                  "Dicamptodontidae",
                                  "Hynobiidae",
                                  "Plethodontidae",
                                  "Proteidae",
                                  "Rhyacotritonidae",
                                  "Salamandridae",
                                  "Sirenidae"),
                       values = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)) +
    
    # add title
    labs(
      x = "Discriminant 1",
      y = "Discriminant 3") +
    
    # custom theme
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.caption = element_text(hjust = 0),
          legend.key = element_rect(fill = NA),
          legend.position = "none", # remove legend
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          # plot.margin = margin(
          #   t = 0.5,
          #   r = 0.5,
          #   b = 0.5,
          #   l = 0.5,
          #   unit = "cm"
          # ),
          aspect.ratio = 1)
  
  dfa_plots[[i]] <- plot1
  dfa_plots[[i+3]] <- plot2
}






##### Phylogenetic PCA ==========================================================================================

### set color palettes

color_lc <- c("#d41159", "#ffc20a", "#1a85ff", "black")



### pca on model fit

# PCA performed on the model fit. That is, on the residuals and covariance of 
# response~predictive not simply mandible. 
# This is different from conventional PCA on centered data
library(mvMORPH)

# loop to compute the analyses
names_fit_list <- names(fit_list_grp_size)
pca_results <- list()

for (i in 1:length(fit_list_grp_size)) {
  model <- fit_list_grp_size[[i]]
  result <- mvgls.pca(model,
                      main = "lc-free PCA")
  
  # store the results
  pca_results[[names_fit_list[i]]] <- result
}

eig_val <- list()
data_pca_eig <- list()

for (i in 1:length(fit_list_grp_size)) {
  # extract pca eigen values
  eig_val[[names_fit_list[i]]] <- as.data.frame(pca_results[[names_fit_list[i]]]$values)
  
  # calculate the percentages
  eig_val[[names_fit_list[i]]]$eig_percent = eig_val[[names_fit_list[i]]] / colSums(eig_val[[names_fit_list[i]]]) * 100
  
  # calculate the cumulative percentage
  eig_val[[names_fit_list[i]]]$eig_cum_percent = cumsum(eig_val[[names_fit_list[i]]]$eig_percent)
}

# extract pca scores
pca_scores <- list()

for (i in 1:length(pca_results)) {
  elements_pca_results <- pca_results[[names_fit_list[i]]]
  
  if ("scores" %in% names(elements_pca_results)) {
    table_scores <- elements_pca_results$scores
    pca_scores[[names_fit_list[i]]] <- table_scores
    
  }
}




### plot
{
  ### create the layout of the plots
  
  mat <- matrix(c(1, 2, 3, 4, 5, 6), # places of the plots (by numeric order)
                nrow = 2,
                ncol = 3)
  layout(mat = mat,
         heights = 1,
         widths = c(1, 1, 1, 1, 1, 1))
  # layout.show(6)
  
  tree <- read.nexus(file = "./data/trees/phylogeny_data.nexus")
  
  color_lc <- setNames(
    c("#d41159", "#ffc20a", "#1a85ff", "black"),
    c("bi", "dd", "pd", "vi"))
  
  
  
  ### 1 - plot pca 1 and 2 all bones
  
  d <- pca_scores[[1]]
  x <- d[,1]
  y <- d[,2]
  
  phylomorphospace(tree = tree,
                   X = d[, 1:2],
                   A = NULL,
                   xlab = "PC 1 - 37%",
                   ylab = "PC 2 - 18.5%",
                   node.size = 0,
                   label = "none", # label of each point
                   lwd = 0.2
  )
  
  points(x = x,
         y = y,
         col = color_lc[lc],
         pch = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)[family],
         cex = 1.5)
  
  
  
  ### 2 - plot pca 1 and 3 all bones
  
  d <- pca_scores[[1]]
  x <- d[,1]
  y <- d[,3]
  
  phylomorphospace(tree = tree,
                   X = d[, c(1, 3)],
                   A = NULL,
                   xlab = "PC 1 - 37%",
                   ylab = "PC 3 - 14.6%",
                   node.size = 0,
                   label = "none", # label of each point
                   lwd = 0.2
  )
  
  points(x = x,
         y = y,
         col = color_lc[lc],
         pch = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)[family],
         cex = 1.5)
  
  
  
  ### 3 - plot pca 1 and 2 dentary
  
  d <- pca_scores[[2]]
  x <- d[,1]
  y <- d[,2]
  
  phylomorphospace(tree = tree,
                   X = d[, c(1, 2)],
                   A = NULL,
                   xlab = "PC 1 - 41.4%",
                   ylab = "PC 2 - 30.7%",
                   node.size = 0,
                   label = "none", # label of each point
                   lwd = 0.2
  )
  
  points(x = x,
         y = y,
         col = color_lc[lc],
         pch = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)[family],
         cex = 1.5)
  
  
  
  ### 4 - plot pca 1 and 3 dentary
  
  d <- pca_scores[[2]]
  x <- d[,1]
  y <- d[,3]
  
  phylomorphospace(tree = tree,
                   X = d[, c(1, 3)],
                   A = NULL,
                   xlab = "PC 1 - 41.4%",
                   ylab = "PC 3 - 6.6%",
                   node.size = 0,
                   label = "none", # label of each point
                   lwd = 0.2
  )
  
  points(x = x,
         y = y,
         col = color_lc[lc],
         pch = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)[family],
         cex = 1.5)
  
  
  
  ### 5 - plot pca 1 and 2 prearticular-angular
  
  d <- pca_scores[[3]]
  x <- d[,1]
  y <- d[,2]
  
  phylomorphospace(tree = tree,
                   X = d[, c(1, 2)],
                   A = NULL,
                   xlab = "PC 1 - 46%",
                   ylab = "PC 2 - 21.8%",
                   node.size = 0,
                   label = "none", # label of each point
                   lwd = 0.2
  )
  
  points(x = x,
         y = y,
         col = color_lc[lc],
         pch = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)[family],
         cex = 1.5)
  
  
  
  ### 6 - plot pca 1 and 3 prearticular-angular
  
  d <- pca_scores[[3]]
  x <- d[,1]
  y <- d[,3]
  
  phylomorphospace(tree = tree,
                   X = d[, c(1, 3)],
                   A = NULL,
                   xlab = "PC 1 - 46%",
                   ylab = "PC 3 - 6.6%",
                   node.size = 0,
                   label = "none", # label of each point
                   lwd = 0.2
  )
  
  points(x = x,
         y = y,
         col = color_lc[lc],
         pch = c(3, 6, 15, 16, 17, 18, 7, 4, 8, 10)[family],
         cex = 1.5)
}







##### visualise allometry =======================================================================================

### extract multivariate scores

# set up color palettes
color_lc <- c("#d41159", "#ffc20a", "#1a85ff", "black")

# plot
names_fit_list_grp_size <- names(fit_list_grp_size)
diag_plot <- list()
par(mfrow = c(3, 3))

library(mvMORPH)
for (i in seq_along(fit_list_grp_size)) {
  
  plot <- plot(fit_list_grp_size[[i]],
               term = "csize",
               residuals = FALSE,
               fitted = TRUE,
               col = color_lc[factor(data_all$lc, levels = c("bi", "dd", "pd", "vi"))],
               pch = 19) 
  abline(h = 0, 
         lty = 2)
  title(main = paste(names(fit_list_grp_size)[i]),
        font.main = 1)
  
  # store the plots
  diag_plot[[names_fit_list_grp_size[i]]] <- plot
}
par(mfrow = c(1, 1))


# extract from plots
load(file = "./data/mandibles/data_lists/data_all.Rdata")
names_fit_list_grp_size <- names(fit_list_grp_size)
mvscores_data <- list()

for (i in 1:length(diag_plot)) {
  mvscores <- diag_plot[[i]]$scores
  
  # assemble data
  mvscores_d <- data.frame(
    mvscores = mvscores,
    lc = data_all$lc,
    lc_fine = data_all$lc_fine,
    habitat = data_all$habitat,
    habitat2 = data_all$habitat2,
    lc_habitat = data_all$lc_habitat,
    csize = data_all$csize
  )
  
  # store the results
  mvscores_data[[names_fit_list_grp_size[i]]] <- mvscores_d
}

names_fit_list_grp_size <- names(fit_list_grp_size)
fitted_grp <- list()
for (i in seq_along(fit_list_grp_size)) {
  res_plot = plot(fit_list_grp_size[[i]],term = "csize",residuals = FALSE,fitted = TRUE) 
  fitted <- fitted(fit_list_grp_size[[i]])%*%res_plot$standBeta
  fitted_grp[[names_fit_list_grp_size[i]]] <- fitted
}

mvscores_data$fit_lc_size_all$fitted <- fitted_grp$fit_lc_size_all
mvscores_data$fit_lc_size_dent$fitted <- fitted_grp$fit_lc_size_dent
mvscores_data$fit_lc_size_artang$fitted <- fitted_grp$fit_lc_size_artang


### plot 

mvscores_plot_lcfine <- list()


library(ggplot2)
for (i in seq_along(mvscores_data)) {
  d <- mvscores_data[[i]]
  
  # set theme
  theme_set(theme_bw())
  
  # plot the data
  plot <- ggplot(d,
                 aes(x = csize, 
                     y = mvscores, 
                     fill = lc_fine,
                     shape = habitat)) +
    geom_point(size = 1.1,
               # shape = 21,
               colour = "black",
               stroke = 0.25,
               alpha = 0.8) +
    geom_line(aes(y = fitted, color = lc, group = lc),
              linewidth = 1) +
    
    # add data on axes
    geom_rug(show.legend = FALSE,
             aes(group = lc_fine, 
                 color = lc_fine)) +
    
    # change legend color and text
    scale_fill_manual(name = "Life cycle",
                      breaks = c("bi", "dd", "pd", "pd1", "pd2", "pd3", "pd4", "vi", "vila", "vipu"),
                      labels = c("bi", "dd", "pd", "pd1", "pd2", "pd3", "pd4", "vi", "vila", "vipu"),
                      values = c("#d41159", "#ffc20a", "#1a85ff", "#85c1e9", "#3498db", "#2874a6", "#1b4f72", "black", "black", "black")) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_color_manual(name = "Life cycle",
                       breaks = c("bi", "dd", "pd", "pd1", "pd2", "pd3", "pd4", "vi", "vila", "vipu"),
                       labels = c("bi", "dd", "pd", "pd1", "pd2", "pd3", "pd4", "vi", "vila", "vipu"),
                       values = c("#d41159", "#ffc20a", "#1a85ff", "#85c1e9", "#3498db", "#2874a6", "#1b4f72", "black", "black", "black")) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_shape_manual(name = "Habitat",
                       breaks = c("A", "SA", "T"),
                       labels = c("Aquatic", "Semi-aquatic", "Terrestrial"),
                       values = c(22, 23, 24)) +
    
    # add title
    labs(x = "Centroid size", y = "Multivariate scores",
         # title = paste("multivariate scores vs centroid size ", names(mvscores_data[i])),
         caption = NULL) + 
    
    # custom theme
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.caption = element_text(hjust = 0),
          legend.key = element_rect(fill = NA),
          # legend.position = "none", # remove legend
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          # plot.margin = margin(
          #   t = 0.5,
          #   r = 0.5,
          #   b = 0.5,
          #   l = 0.5,
          #   unit = "cm"
          # ),
          aspect.ratio = 1)
  
  # store the plots
  mvscores_plot_lcfine[[i]] <- plot
}







##### disparity ==============================================================================================

load(file = "./data/mandibles/data_lists/data_all.Rdata")
load(file = "./data/mandibles/data_lists/data_dent.Rdata")
load(file = "./data/mandibles/data_lists/data_artang.Rdata")
load(file = "./data/mandibles/models_fits/fit_list_grp_size.Rdata")

names_fit_list <- names(fit_list_grp_size)
shape_data <- list()

for (i in 1:length(fit_list_grp_size)) {
  model <- fit_list_grp_size[[i]]
  
  residuals_val <- residuals(model,
                             type = "normalized")
  
  fitted_val <- fitted(model)
  
  new_shape_data <- residuals_val + fitted_val
  rownames(new_shape_data) <- sp_names
  
  # store the results
  shape_data[[names_fit_list[i]]] <- new_shape_data
}



### making groups for analyses

lc_habitat_df <- as.data.frame(lc_habitat)
rownames(lc_habitat_df) <- sp_names



### dispRity 

library(dispRity)

names_fit_list <- names(fit_list_grp_size)
disparity_results <- list()

for (i in 1:length(shape_data)) {
  
  d_shape <- shape_data[[i]]
  
  # disparity for life cycles
  result <- dispRity.per.group(data = d_shape,
                               group = lc_habitat_df,
                               metric = c(sum, variances))
  
  # store the results
  disparity_results[[names_fit_list[i]]] <- result
}



### post hoc test 

test_disp_lc_habitat <- list()

library(dispRity)

for (i in 1:length(disparity_results)) {
  # test disparity
  test_disp <- test.dispRity(disparity_results[[i]],
                             test = wilcox.test,
                             comparisons = "pairwise",
                             correction = "bonferroni")
  
  test_disp_df <- as.data.frame(test_disp)
  
  # post hoc for life cycles
  test_disp_lc_habitat[[names_fit_list[i]]] <- test_disp_df
}



### plot

disp_obs <- get.disparity(disparity_results$fit_lc_size_artang)

# prepare data

observed_disparity_lc_habitat <- list()

for (i in 1:length(disparity_results)) {

  # extract disparity scores
  disp_obs <- get.disparity(disparity_results[[i]],
                            observed = TRUE)
  
  disp_obs_df <- as.data.frame(disp_obs)
  
  observed_disparity_lc_habitat[[names_fit_list[i]]] <- disp_obs_df
  
  
  # assemble into 1 data frame
  disp_lc <- do.call(rbind, observed_disparity_lc_habitat)
  
  # add a column with bone identification
  disp_lc$bone <- substr(rownames(disp_lc),
                         13,
                         nchar(rownames(disp_lc)))
  
  # renames columns
  colnames(disp_lc) <- c("bi_A",
                         "bi_SA",
                         "bi_T",
                         "dd_T",
                         "pd_A",
                         "vi_T",
                         "bone")
  
  # remove rownames
  rownames(disp_lc) <- NULL
  
  # transpose data frame
  library(reshape2)
  data_disp_lc <- melt(disp_lc)
  
  # renames columns
  colnames(data_disp_lc) <- c("bone",
                              "lc",
                              "observed_disparity")
  
  # arrangge data frame
  library(dplyr)
  data_disp_lc <- data_disp_lc  %>% 
    arrange(match(bone,
                  c("all",
                    "dent",
                    "artang")))  # order
  
  data_disp_lc$bone <- factor(data_disp_lc$bone,
                              levels = c("all",
                                         "dent",
                                         "artang")) # lock in factor level order
}
list_data_disp <- data_disp_lc



# modify y axis for disparity plots
list_data_disp_fig <- list_data_disp

list_data_disp_fig$observed_disparity <- list_data_disp_fig$observed_disparity * 100

# ensure that ggplot will plot data in the correct order
library(dplyr)
list_data_disp_fig$bone <- factor(list_data_disp_fig$bone, 
                                  levels = c("all", "dent", "artang"))
list_data_disp_fig$lc <- factor(list_data_disp_fig$lc, 
                                levels = c("bi_A", "bi_SA", "bi_T", "dd_T", "pd_A", "vi_T"))

# plot disparity with new values on y axis
disp_plots_fig <- list()

d <- list_data_disp_fig
# set theme
theme_set(theme_bw())
plot <- ggplot(d,
               aes(x = bone,
                   y = observed_disparity,
                   fill = lc,
                   shape = lc)) +
  geom_point(size = 2.5,
             # shape = 21,
             colour = "black",
             stroke = 0.25,
             alpha = 0.8) +
  
  scale_fill_manual(name = "Life cycle",
                    breaks = c("bi_A", "bi_SA", "bi_T", "dd_T", "pd_A", "vi_T"),
                    labels = c("biphasic", "biphasic", "biphasic", "direct development", "paedomorphic", "viviparous"),
                    values = c("#d41159", "#d41159", "#d41159", "#ffc20a", "#1a85ff", "black")) +
  scale_shape_manual(name = "Habitat",
                     breaks = c("bi_A", "bi_SA", "bi_T", "dd_T", "pd_A", "vi_T"),
                     labels = c("Aquatic", "Semi-aquatic", "Terrestrial", "Terrestrial", "Aquatic", "Terrestrial"),
                     values = c(22, 21, 24, 24, 22, 24)) +
  
  # add title
  labs(x = NULL, 
       y = "Disparity",
       # title = "Disparity in mandible shape in Caudata \nwith regard to life cycle",
       caption = NULL) +
  
  # change labels on x axis
  scale_x_discrete(labels = c("All", 
                              "Dent", 
                              "Artang")) +
  
  # custom theme
  theme(plot.title = element_text(size = 9, face = "bold"),
        plot.caption = element_text(hjust = 0),
        legend.key = element_rect(fill = NA),
        # legend.position = "none", # remove legend
        # axis.title = element_text(size = 8),
        # axis.text = element_text(size = 8),
        # axis.title.x = element_blank(), # remove axis names
        # axis.title.y = element_blank(), # remove axis names
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        # plot.margin = margin(
        #   t = 0.5,
        #   r = 0.5,
        #   b = 0.5,
        #   l = 0.5,
        #   unit = "cm"
        # ),
        aspect.ratio = 1)

disp_plots_fig <- plot







##### rate of evolution ========================================================================================

tree = phy

data_lc_habitat <- data_all$lc_habitat

nsim = 1000 # number of stochastic mapping

tree$branch.length <- tree$branch.length/max(nodeHeights(tree)) # standardise tree branches
my_trees <- make.simmap(tree, 
                        data_lc, 
                        model = "ARD", 
                        nsim = nsim)
distrib <- summary(my_trees,
                   plot = TRUE)

# running models
library(mvMORPH)

all_rates_all <- sapply(1:nsim, function(x) {
  # model fit
  fit_trees <- mvgls(Y ~ lc * csize,
                     data = data_lc_habitat,
                     tree = my_trees[[x]],
                     model = "BMM",
                     method = "PL",
                     error = TRUE)
  fit_trees$param # estimate rates
})

all_rates_dent <- sapply(1:nsim, function(x) {
  # model fit
  fit_trees <- mvgls(Y ~ lc * csize,
                     data = data_lc_habitat,
                     tree = my_trees[[x]],
                     model = "BMM",
                     method = "PL",
                     error = TRUE)
  fit_trees$param # estimate rates
})

all_rates_artang <- sapply(1:nsim, function(x) {
  # model fit
  fit_trees <- mvgls(Y ~ lc * csize,
                     data = data_lc_habitat,
                     tree = my_trees[[x]],
                     model = "BMM",
                     method = "PL",
                     error = TRUE)
  fit_trees$param # estimate rates
})





### post-hoc tests

# calculate rowmeans
rate_all <-colMeans(all_rates_all)
rate_dent <-colMeans(all_rates_dent)
rate_artang <-colMeans(all_rates_artang)



### assess differences of means

data <- all_rates_all
data <- all_rates_dent
data <- all_rates_artang

# initialise matrice of p values
p_values <- matrix(nrow = ncol(data), ncol = ncol(data))
colnames(p_values) <- rownames(p_values) <- colnames(data)

# test for differences between groups
total_comparisons <- ncol(data) * (ncol(data) - 1) / 2
for (i in 1:ncol(data)) {
  for (j in 1:ncol(data)) {
    if (i == j) {
      p_values[i, j] <- NA  # don't consider comparison of the same variable
    } else {
      test_result <- wilcox.test(data[, i], data[, j])
      # apply Bonferroni correction
      p_values[i, j] <- min(test_result$p.value * total_comparisons, 1)
    }
  }
}

rates_all_pvalue <- p_values
rates_dent_pvalue <- p_values
rates_artang_pvalue <- p_values


### plot 

# assemble data frames
rate_sim <- rbind(rate_all,
                  rate_dent, 
                  rate_artang)

# for some reason, you need to transform it again into dataframe
rate_sim <- as.data.frame(rate_sim)

# add a column bone
rate_sim$bone <- c("all",
                   "dentary",
                   "articular-angular")

# transform data
library(reshape2)
rate_sim_dat <- melt(rate_sim)

# rename columns of new data frame
colnames(rate_sim_dat) <- c("bone",
                            "life_cycle",
                            "mean_rate")

# arrange data 
library(dplyr)
rate_sim_dat <- rate_sim_dat  %>% 
  arrange(match(bone, 
                c("all", "dentary",  "articular-angular")))  # order
rate_sim_dat$bone <- factor(rate_sim_dat$bone,
                            levels = c("all", "dentary",  "articular-angular")) # lock in factor level order


# modify y axis for rate of evolution plots
rate_sim_dat_fig <- rate_sim_dat
rate_sim_dat_fig$mean_rate <- rate_sim_dat_fig$mean_rate * 100000

# ensure that ggplot will plot data in the correct order
library(dplyr)
rate_sim_dat_fig$bone <- factor(rate_sim_dat_fig$bone, 
                                levels = c("all", "dentary", "articular-angular"))
rate_sim_dat_fig$lc <- factor(rate_sim_dat_fig$lc, 
                              levels = c("bi_A", "bi_SA", "bi_T", "dd_T", "pd_A", "vi_T"))


# plot rate of evolution with new values on y axis
library(ggplot2)
# set them
theme_set(theme_bw())
# plot the data
plot_roe_fig <- ggplot(rate_sim_dat_fig,
                       aes(x = bone,
                           y = mean_rate,
                           fill = life_cycle, 
                           shape = life_cycle)) +
  geom_point(size = 2.5,
             colour = "black",
             stroke = 0.25,
             alpha = 0.8) +
  # change legend color and text
  scale_fill_manual(name = "Life cycle",
                    breaks = c("bi_A", "bi_SA", "bi_T", "dd_T", "pd_A", "vi_T"),
                    labels = c("biphasic", "biphasic", "biphasic", "direct development", "paedomorphic", "viviparous"),
                    values = c("#d41159", "#d41159", "#d41159", "#ffc20a", "#1a85ff", "black")) +
  scale_shape_manual(name = "Habitat",
                     breaks = c("bi_A", "bi_SA", "bi_T", "dd_T", "pd_A", "vi_T"),
                     labels = c("Aquatic", "Semi-aquatic", "Terrestrial", "Terrestrial", "Aquatic", "Terrestrial"),
                     values = c(22, 21, 24, 24, 22, 24)) +
  # add title
  labs(x = NULL,
       y = "Rate of evolution",
       caption = NULL) +
  # change x labels
  scale_x_discrete(labels = c("All", 
                              "Dent", 
                              "Artang")) +
  scale_y_continuous(labels = function(mean_rate) sprintf("%.1f", mean_rate)) +
  
  # custom theme
  theme(plot.title = element_text(size = 9, face = "bold"),
        plot.caption = element_text(hjust = 0),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        aspect.ratio = 1)








##### integration and modularity ===============================================================================

### subset data with lc and habitat

# identify species for each life cycle
sp_bi_A <- subset(data_all$lc_habitat,
                  data_all$lc_habitat == "bi_A")
spbi_A <- names(sp_bi_A)

sp_bi_SA <- subset(data_all$lc_habitat,
                   data_all$lc_habitat == "bi_SA")
spbi_SA <- names(sp_bi_SA)

sp_bi_T <- subset(data_all$lc_habitat,
                  data_all$lc_habitat == "bi_T")
spbi_T <- names(sp_bi_T)

sp_dd_T <- subset(data_all$lc_habitat,
                  data_all$lc_habitat == "dd_T")
spdd_T <- names(sp_dd_T)

sp_pd_A <- subset(data_all$lc_habitat,
                  data_all$lc_habitat == "pd_A")
sppd_A <- names(sp_pd_A)

sp_vi_T <- subset(data_all$lc_habitat,
                  data_all$lc_habitat == "vi_T")
spvi_T <- names(sp_vi_T)

# identify landmarks for each bone
dentary <- c(1:4, 7, 10:44, 61:70)
artang <- c(5:6, 45:60)
coronoid <- c(8:9, 71:80)

# subset landmarks for eahc life cycle and remove coronoid
load(file = "./data/mandibles/landmarks/Y.Rdata")

Y_bi_A_no_cor <- Y[-coronoid, , spbi_A]
Y_bi_SA_no_cor <- Y[-coronoid, , spbi_SA]
Y_bi_T_no_cor <- Y[-coronoid, , spbi_T]
Y_dd_T_no_cor <- Y[-coronoid, , spdd_T]
Y_pd_A_no_cor <- Y[-coronoid, , sppd_A]
Y_vi_T_no_cor <- Y[-coronoid, , spvi_T]

# create partition for each bone 
partition <- c("d", "d", "d", "d",
               "a", "a",
               "d", 
               "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d", "d",
               "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", 
               "d", "d", "d", "d", "d", "d", "d", "d", "d", "d")



### create trees that matches the subset data 

# get phylogeny
library(geiger)
library(ape)
phy <- read.nexus(file="./data/trees/MCC_tree_graft_pd.nexus")

# 2 d arrays
library(geomorph)
Y_bi_A_no_cor_2d <- two.d.array(Y_bi_A_no_cor)
Y_bi_SA_no_cor_2d <- two.d.array(Y_bi_SA_no_cor)
Y_bi_T_no_cor_2d <- two.d.array(Y_bi_T_no_cor)
Y_dd_T_no_cor_2d <- two.d.array(Y_dd_T_no_cor)
Y_pd_A_no_cor_2d <- two.d.array(Y_pd_A_no_cor)
Y_vi_T_no_cor_2d <- two.d.array(Y_vi_T_no_cor)

# match the phylogeny with dataset
td_bi_A_no_cor_2d <- treedata(phy, Y_bi_A_no_cor_2d) 
td_bi_SA_no_cor_2d <- treedata(phy, Y_bi_SA_no_cor_2d) 
td_bi_T_no_cor_2d <- treedata(phy, Y_bi_T_no_cor_2d) 
td_dd_T_no_cor_2d <- treedata(phy, Y_dd_T_no_cor_2d) 
td_pd_A_no_cor_2d <- treedata(phy, Y_pd_A_no_cor_2d) 
td_vi_T_no_cor_2d <- treedata(phy, Y_vi_T_no_cor_2d) 

# extract phylogy
phy_bi_A_no_cor = td_bi_A_no_cor_2d$phy 
phy_bi_SA_no_cor = td_bi_SA_no_cor_2d$phy 
phy_bi_T_no_cor = td_bi_T_no_cor_2d$phy 
phy_dd_T_no_cor = td_dd_T_no_cor_2d$phy 
phy_pd_A_no_cor = td_pd_A_no_cor_2d$phy
phy_vi_T_no_cor = td_vi_T_no_cor_2d$phy



### Integration

library(geomorph)

# using data in 3d array
{
  int_bi_A <- phylo.integration(A = Y_bi_A_no_cor,
                                partition.gp = partition,
                                phy = phy_bi_A_no_cor,
                                iter = 999,
                                seed = "random",
                                print.progress = TRUE)
  
  int_bi_SA <- phylo.integration(A = Y_bi_SA_no_cor,
                                 partition.gp = partition,
                                 phy = phy_bi_SA_no_cor,
                                 iter = 999,
                                 seed = "random",
                                 print.progress = TRUE)
  
  int_bi_T <- phylo.integration(A = Y_bi_T_no_cor,
                                partition.gp = partition,
                                phy = phy_bi_T_no_cor,
                                iter = 999,
                                seed = "random",
                                print.progress = TRUE)
  
  int_dd_T <- phylo.integration(A = Y_dd_T_no_cor,
                                partition.gp = partition,
                                phy = phy_dd_T_no_cor,
                                iter = 999,
                                seed = "random",
                                print.progress = TRUE)
  
  int_pd_A <- phylo.integration(A = Y_pd_A_no_cor,
                                partition.gp = partition,
                                phy = phy_pd_A_no_cor,
                                iter = 999,
                                seed = "random",
                                print.progress = TRUE)
  
  int_vi_T <- phylo.integration(A = Y_vi_T_no_cor,
                                partition.gp = partition,
                                phy = phy_vi_T_no_cor,
                                iter = 999,
                                seed = "random",
                                print.progress = TRUE)
}

# save results
table_integration <- data.frame(
  life_cycle = c("bi_a",
                 "bi_sa",
                 "bi_t",
                 "dd",
                 "pd",
                 "vi"),
  r_pls = c(int_bi_A$r.pls,
            int_bi_SA$r.pls,
            int_bi_T$r.pls,
            int_dd_T$r.pls,
            int_pd_A$r.pls,
            int_vi_T$r.pls),
  Effect_size = c(int_bi_A$Z,
                  int_bi_SA$Z,
                  int_bi_T$Z,
                  int_dd_T$Z,
                  int_pd_A$Z,
                  int_vi_T$Z),
  P_value = c(int_bi_A$P.value,
              int_bi_SA$P.value,
              int_bi_T$P.value,
              int_dd_T$P.value,
              int_pd_A$P.value,
              int_vi_T$P.value)
)



### plot

table_integration <- read.csv2(file = "./results/integration_modularity/table_integration.csv")
table_integration_plot <- table_integration
table_integration_plot$habitat <- c("a", "sa", "t", "t", "a", "t")
table_integration_plot$all_sp <- c("all_sp", "all_sp", "all_sp", "all_sp", "all_sp", "all_sp")

# plot rate of evolution with new values on y axis
library(ggplot2)
# set them
theme_set(theme_bw())
# plot the data
plot_int_fig <- ggplot(table_integration_plot,
                       aes(x = all_sp,
                           y = r_pls,
                           fill = life_cycle,
                           shape = habitat)) +
  
  geom_point(size = 2.5,
             colour = "black",
             stroke = 0.25,
             alpha = 0.8) +
  
  scale_fill_manual(name = "Life cycle",
                    breaks = c("bi_a", "bi_sa", "bi_t", "dd", "pd", "vi"),
                    labels = c("Biphasic aquatic", "Biphasic semi-aquatic", "Biphasic terrestrial", "Direct development", "Paedomorphic", "Viviparous"),
                    values = c("#d41159", "#d41159", "#d41159", "#ffc20a", "#1a85ff", "black")) +
  
  scale_shape_manual(name = "Habitat",
                     breaks = c("a", "sa", "t", "t", "a", "t"),
                     labels = c("Aquatic", "Semi-aquatic", "Terrestrial", "Terrestrial", "Aquatic", "Terrestrial"),
                     values = c(22, 21, 24, 24, 22, 24)) +
  
  labs(x = NULL,
       y = "Integration (R-PLS)",
       caption = NULL) +
  
  # custom theme
  theme(plot.title = element_text(size = 9, face = "bold"),
        plot.caption = element_text(hjust = 0),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        aspect.ratio = 1)






### Modularity

library(geomorph)
{
  modul_bi_A <- phylo.modularity(A = Y_bi_A_no_cor,
                                 partition.gp = partition,
                                 phy = phy_bi_A_no_cor,
                                 CI = TRUE,
                                 iter = 999,
                                 seed = "random",
                                 print.progress = TRUE)
  
  modul_bi_SA <- phylo.modularity(A = Y_bi_SA_no_cor,
                                  partition.gp = partition,
                                  phy = phy_bi_SA_no_cor,
                                  CI = TRUE,
                                  iter = 999,
                                  seed = "random",
                                  print.progress = TRUE)
  
  modul_bi_T <- phylo.modularity(A = Y_bi_T_no_cor,
                                 partition.gp = partition,
                                 phy = phy_bi_T_no_cor,
                                 CI = TRUE,
                                 iter = 999,
                                 seed = "random",
                                 print.progress = TRUE)
  
  modul_dd_T <- phylo.modularity(A = Y_dd_T_no_cor,
                                 partition.gp = partition,
                                 phy = phy_dd_T_no_cor,
                                 CI = TRUE,
                                 iter = 999,
                                 seed = "random",
                                 print.progress = TRUE)
  
  modul_pd_A <- phylo.modularity(A = Y_pd_A_no_cor,
                                 partition.gp = partition,
                                 phy = phy_pd_A_no_cor,
                                 CI = TRUE,
                                 iter = 999,
                                 seed = "random",
                                 print.progress = TRUE)
  
  modul_vi_T <- phylo.modularity(A = Y_vi_T_no_cor,
                                 partition.gp = partition,
                                 phy = phy_vi_T_no_cor,
                                 CI = TRUE,
                                 iter = 999,
                                 seed = "random",
                                 print.progress = TRUE)
}

# save results
table_modularity <- data.frame(
  row.names =  c("Biphasic aquatic",
                 "Biphasic semi-aquatic",
                 "Biphasic terrestrial",
                 "Direct development terrestrial",
                 "Paedomorphic aquatic",
                 "Viviparous terrestrial"),
  life_cycle = c("bi_a",
                 "bi_sa",
                 "bi_t",
                 "dd",
                 "pd",
                 "vi"),
  cr = c(modul_bi_A$CR,
         modul_bi_SA$CR,
         modul_bi_T$CR,
         modul_dd_T$CR,
         modul_pd_A$CR,
         modul_vi_T$CR),
  ci_2.5 = c(modul_bi_A$CInterval[[1]],
             modul_bi_SA$CInterval[[1]],
             modul_bi_T$CInterval[[1]],
             modul_dd_T$CInterval[[1]],
             modul_pd_A$CInterval[[1]],
             modul_vi_T$CInterval[[1]]),
  ci_97.5 = c(modul_bi_A$CInterval[[2]],
              modul_bi_SA$CInterval[[2]],
              modul_bi_T$CInterval[[2]],
              modul_dd_T$CInterval[[2]],
              modul_pd_A$CInterval[[2]],
              modul_vi_T$CInterval[[2]]),
  effect_size = c(modul_bi_A$Z,
                  modul_bi_SA$Z,
                  modul_bi_T$Z,
                  modul_dd_T$Z,
                  modul_pd_A$Z,
                  modul_vi_T$Z),
  p_value = c(modul_bi_A$P.value,
              modul_bi_SA$P.value,
              modul_bi_T$P.value,
              modul_dd_T$P.value,
              modul_pd_A$P.value,
              modul_vi_T$P.value)
)


### plot results


table_modularity_fig <- table_modularity
table_modularity_fig$species <- c("all", "all", "all", "all", "all", "all")

# save table
write.table(table_modularity_fig,
            file = paste0("./results/table_modularity_fig.txt"),
            sep = ";")

library(ggplot2)
# set them
theme_set(theme_bw())
# plot the data
plot_cr_all_bi_fig <- ggplot(table_modularity_fig,
                             aes(x = species,
                                 y = cr,
                                 fill = life_cycle,
                                 shape = life_cycle)) +
  geom_point(size = 2.5,
             colour = "black",
             stroke = 0.25,
             alpha = 0.8) +
  
  # change legend color and text
  scale_fill_manual(name = "Life cycle",
                    breaks = c("bi_a", "bi_sa", "bi_t", "dd", "pd", "vi"),
                    labels = c("biphasic", "biphasic", "biphasic", "direct development", "paedomorphic", "viviparous"),
                    values = c("#d41159", "#d41159", "#d41159", "#ffc20a", "#1a85ff", "black")) +
  scale_shape_manual(name = "Habitat",
                     breaks = c("bi_a", "bi_sa", "bi_t", "dd", "pd", "vi"),
                     labels = c("Aquatic", "Semi-aquatic", "Terrestrial", "Terrestrial", "Aquatic", "Terrestrial"),
                     values = c(22, 21, 24, 24, 22, 24)) +
  
  # add title
  labs(x = NULL,
       y = "Modularity (CR)",
       caption = NULL) +

  # custom theme
  theme(plot.title = element_text(size = 9, face = "bold"),
        plot.caption = element_text(hjust = 0),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        aspect.ratio = 1)






##### ancestral state reconstruction ==========================================================================

### phylogeny

library(ape)
caudata_tree <- read.nexus(file = "./data/trees/phylogeny_data.nexus")
caudata_phy <- ladderize(caudata_tree)



### species data

caudata_data <- read.csv2(file = "./data/mandibles/ident_tables/ident_ordered.csv")

# extract characters as vectors
lc_4 <- setNames(caudata_data$lc,
                 caudata_data$sp)



### fit Mk model with phytools

tree <- caudata_phy
x <- lc_4



## Equal rates (for permitted transitions) model ER

# multiple iterations
library(phytools)
library(foreach)
library(doParallel)
niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl = mc)

all_er_fits <- foreach(i = 1:niter) %dopar% {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error")) {
    obj <- try(phytools::fitMk(
      tree = tree,
      x = x,
      model = "ER",
      pi = "equal",
      logscale = sample(c(FALSE, TRUE), 1),
      opt.method = sample(c("nlminb", "optim"), 1),
      rand.start = TRUE
    ))
  }
  obj
}
stopCluster(mc) ## stop cluster

# best fitting model
lnL <- sapply(all_er_fits, logLik)
fit_best_er <- all_er_fits[[which.max(lnL)]]



## symetric backward and foward rate (for permitted transitions) model SYM

# multiple iterations
library(phytools)
library(foreach)
library(doParallel)
niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl = mc)

all_sym_fits <- foreach(i = 1:niter) %dopar% {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error")) {
    obj <- try(phytools::fitMk(
      tree = tree,
      x = x,
      model = "SYM",
      pi = "equal",
      logscale = sample(c(FALSE, TRUE), 1),
      opt.method = sample(c("nlminb", "optim"), 1),
      rand.start = TRUE
    ))
  }
  obj
}
stopCluster(mc) ## stop cluster

# best fitting model
lnL <- sapply(all_sym_fits, logLik)
fit_best_sym <- all_sym_fits[[which.max(lnL)]]



## different rates (for permitted transitions) model ARD

# multiple iterations
library(phytools)
library(foreach)
library(doParallel)
niter <- 1000 # set iterations
# set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl = mc)

all_ard_fits <- foreach(i = 1:niter) %dopar% {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error")) {
    obj <- try(phytools::fitMk(
      tree = tree,
      x = x,
      model = "ARD",
      pi = "equal",
      logscale = sample(c(FALSE, TRUE), 1),
      opt.method = sample(c("nlminb", "optim"), 1),
      rand.start = TRUE
    ))
  }
  obj
}
stopCluster(mc) ## stop cluster

# best fitting model
lnL <- sapply(all_ard_fits, logLik)
fit_best_ard <- all_ard_fits[[which.max(lnL)]]





### fitting custom models

tree <- caudata_tree
x <- lc_4

# all transitions possible except for: vi -> bi / dd / pd
custom_model_vi <- matrix(c(
  0,1,2,3,
  4,0,5,6,
  7,8,0,9,
  0,0,0,0), 4, 4,
  byrow = TRUE,
  dimnames = list(c("bi", "dd", "pd", "vi"), 
                  c("bi", "dd", "pd", "vi")))
custom_model_vi

# all transitions possible except for: vi -> bi / dd / pd and pd -> bi / dd / vi
custom_model_vi_pd <- matrix(c(
  0,1,2,3,
  4,0,5,6,
  0,0,0,0,
  0,0,0,0), 4, 4,
  byrow = TRUE,
  dimnames = list(c("bi", "dd", "pd", "vi"), 
                  c("bi", "dd", "pd", "vi")))
custom_model_vi_pd

# all transitions possible except for: pd -> bi / dd / vi
custom_model_pd <- matrix(c(
  0,1,2,3,
  4,0,5,6,
  0,0,0,0,
  7,8,0,9), 4, 4,
  byrow = TRUE,
  dimnames = list(c("bi", "dd", "pd", "vi"), 
                  c("bi", "dd", "pd", "vi")))
custom_model_pd

# run multiple optimization iterations to try to ensure that we converge to the true Maximum Likelihood solution
library(phytools)
library(foreach)
library(doParallel)
niter <- 1000 ## set iterations
## set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl = mc)

all_custom_vi_fits <- foreach(i = 1:niter) %dopar% {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error")) {
    obj <- try(phytools::fitMk(
      tree,
      x,
      model = custom_model_vi,
      pi = "equal",
      logscale = sample(c(FALSE, TRUE), 1),
      opt.method = sample(c("nlminb", "optim"), 1),
      rand.start = TRUE
    ))
  }
  obj
}
stopCluster(mc) ## stop cluster

# best fitting model
lnL <- sapply(all_custom_vi_fits, logLik)
fit_best_custom_vi <- all_custom_vi_fits[[which.max(lnL)]]


library(phytools)
library(foreach)
library(doParallel)
niter <- 1000 ## set iterations
## set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl = mc)

all_custom_vi_pd_fits <- foreach(i = 1:niter) %dopar% {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error")) {
    obj <- try(phytools::fitMk(
      tree,
      x,
      model = custom_model_vi_pd,
      pi = "equal",
      logscale = sample(c(FALSE, TRUE), 1),
      opt.method = sample(c("nlminb", "optim"), 1),
      rand.start = TRUE
    ))
  }
  obj
}
stopCluster(mc) ## stop cluster

# best fitting model
lnL <- sapply(all_custom_vi_pd_fits, logLik)
fit_best_custom_vi_pd <- all_custom_vi_pd_fits[[which.max(lnL)]]


library(phytools)
library(foreach)
library(doParallel)
niter <- 1000 ## set iterations
## set ncores and open cluster
ncores <- min(niter, detectCores() - 1)
mc <- makeCluster(ncores, type = "PSOCK")
registerDoParallel(cl = mc)

all_custom_pd_fits <- foreach(i = 1:niter) %dopar% {
  obj <- list()
  class(obj) <- "try-error"
  while (inherits(obj, "try-error")) {
    obj <- try(phytools::fitMk(
      tree,
      x,
      model = custom_model_pd,
      pi = "equal",
      logscale = sample(c(FALSE, TRUE), 1),
      opt.method = sample(c("nlminb", "optim"), 1),
      rand.start = TRUE
    ))
  }
  obj
}
stopCluster(mc) ## stop cluster

# best fitting model
lnL <- sapply(all_custom_pd_fits, logLik)
fit_best_custom_pd <- all_custom_pd_fits[[which.max(lnL)]]



### compare the fitted models 

# results as data frame
rslt_anova <- anova(fit_best_er,
                    fit_best_sym,
                    fit_best_ard,
                    fit_best_custom_pd,
                    fit_best_custom_vi,
                    fit_best_custom_vi_pd)

rslt_anc_models <- data.frame(
  model = rownames(rslt_anova),
  logLik = rslt_anova[,1],
  df = rslt_anova$d.f.,
  AIC = rslt_anova$AIC,
  weight = rslt_anova$weight)
library(dplyr)
rslt_anc_models <- rslt_anc_models %>% mutate(Delta_AIC = AIC - min(AIC))
rslt_anc_models <- rslt_anc_models[order(rslt_anc_models$Delta_AIC), ]
rslt_anc_models



### ancestral states reconstruction 

library(phytools)
fit_sym_anc <- ancr(fit_best_sym,
                    type = "joint")



### plot model rates 

cols = c("#d41159", "#ffc20a", "#1a85ff", "black")

pdf(file = "./results/ancestral_state_reconstruction/figure_caudata_asr_all_model_rates.pdf",
    height = 7,
    width = 8)

{
  mat <- matrix(c(1,4,2,5,3,6), 
                nrow = 2,
                ncol = 3)
  
  layout(mat)
  # layout.show(n = 6)
  
  library(plotrix)
  
  # er
  fit_mod <- fit_best_er
  p <- plot(
    fit_mod,
    width = TRUE,
    color = FALSE,
    text = TRUE,
    show.zeros = TRUE,
    max.lwd = 6,
    cex.traits = 1,
    spacer = 0.15,
    offset = 0.05
  )
  invisible(mapply(plotrix::draw.circle, 
                   x = p$x, 
                   y = p$y, 
                   col = cols, 
                   MoreArgs = list(radius = strheight("0"), border = "transparent")))
  text(p$x,
       p$y,
       p$states,
       col = c("white", "white", "white", "white"))
  mtext("ER",
        line = 0,
        adj = 0.5)
  
  # sym
  fit_mod <- fit_best_sym
  p <- plot(
    fit_mod,
    width = TRUE,
    color = FALSE,
    text = TRUE,
    show.zeros = TRUE,
    max.lwd = 6,
    cex.traits = 1, 
    spacer = 0.15, 
    offset = 0.05
  )
  invisible(mapply(plotrix::draw.circle, 
                   x = p$x, 
                   y = p$y, 
                   col = cols, 
                   MoreArgs = list(radius = strheight("0"), border = "transparent")))
  text(p$x,
       p$y,
       p$states,
       col = c("white", "white", "white", "white"))
  mtext("SYM",
        line = 0,
        adj = 0.5)
  
  # ard
  fit_mod <- fit_best_ard
  p <- plot(
    fit_mod,
    width = TRUE,
    color = FALSE,
    text = TRUE,
    show.zeros = TRUE,
    max.lwd = 6,
    cex.traits = 1, 
    spacer = 0.15, 
    offset = 0.05
  )
  invisible(mapply(plotrix::draw.circle, 
                   x = p$x, 
                   y = p$y, 
                   col = cols, 
                   MoreArgs = list(radius = strheight("0"), border = "transparent")))
  text(p$x,
       p$y,
       p$states,
       col = c("white", "white", "white", "white"))
  mtext("ARD",
        line = 0,
        adj = 0.5)
  
  # custom vi
  fit_mod <- fit_best_custom_vi
  p <- plot(
    fit_mod, 
    width = TRUE, 
    color = FALSE, 
    text = TRUE, 
    show.zeros = TRUE, 
    max.lwd = 6, 
    cex.traits = 1, 
    spacer = 0.15, 
    offset = 0.05
  )
  invisible(mapply(plotrix::draw.circle, 
                   x = p$x, 
                   y = p$y, 
                   col = cols, 
                   MoreArgs = list(radius = strheight("0"), border = "transparent")))
  text(p$x,
       p$y,
       p$states,
       col = c("white", "white", "white", "white"))
  mtext("Viviparous dead end",
        line = 0,
        adj = 0.5)
  
  # custom pd
  fit_mod <- fit_best_custom_pd
  p <- plot(
    fit_mod, 
    width = TRUE, 
    color = FALSE, 
    text = TRUE, 
    show.zeros = TRUE, 
    max.lwd = 6, 
    cex.traits = 1, 
    spacer = 0.15, 
    offset = 0.05
  )
  invisible(mapply(plotrix::draw.circle, 
                   x = p$x, 
                   y = p$y, 
                   col = cols, 
                   MoreArgs = list(radius = strheight("0"), border = "transparent")))
  text(p$x,
       p$y,
       p$states,
       col = c("white", "white", "white", "white"))
  mtext("Paedormorph dead end",
        line = 0,
        adj = 0.5)
  
  # custom vi pd
  fit_mod <- fit_best_custom_vi_pd
  p <- plot(
    fit_mod, 
    width = TRUE, 
    color = FALSE, 
    text = TRUE, 
    show.zeros = TRUE, 
    max.lwd = 6, 
    cex.traits = 1, 
    spacer = 0.15, 
    offset = 0.05
  )
  invisible(mapply(plotrix::draw.circle, 
                   x = p$x, 
                   y = p$y, 
                   col = cols, 
                   MoreArgs = list(radius = strheight("0"), border = "transparent")))
  text(p$x,
       p$y,
       p$states,
       col = c("white", "white", "white", "white"))
  mtext("Viviparous and Paedomorph dead end",
        line = 0,
        adj = 0.5)
}
dev.off()



### plot ancestral reconstruction

load(file = "./results/ancestral_state_reconstruction/fit_sym_anc.Rdata")

cols <- c("#d41159", "#ffc20a", "#1a85ff", "black")

fit_mod_anc <- fit_sym_anc

cex = 0.1

plot(fit_sym_anc,
     type = "arc",
     fsize = 0.1,
     direction = "rightwards",
     args.nodelabels = list(cex = cex, piecol = cols),
     arc_height = 0.1,
     mar = c(0.1, 0.1, 1.1, 0.1))
mtext("Model marginal ancestral states",
      line = 0,
      adj = 0)



### plot lineages through time

# extract plot data
library(phytools)

d_ave_sim_sym_ltt <- read.csv2(file = "./results/ancestral_state_reconstruction/d_ave_sim_sym_ltt.csv")
{
  # plot framework
  par(mar = c(4, 4, 4, 4))
  plot(1,
       type = 'n',
       ylim = c(round(min(d_ave_sim_sym_ltt$dd), 1), round(max(d_ave_sim_sym_ltt$dd), 1) + 10),
       xlim = c(round(min(d_ave_sim_sym_ltt$time), 1) - 6, round(max(d_ave_sim_sym_ltt$time), 1)), 
       axes = FALSE, 
       ylab = '',
       xlab = '',
       main = '',
       bty = 'n',
       las = 1,
       cex.axis = 0.5,
       cex.lab = 0.5)
  
  # lines
  lines(d_ave_sim_sym_ltt$time,
        d_ave_sim_sym_ltt$bi,
        col = "#d41159", 
        lwd = 2)
  lines(d_ave_sim_sym_ltt$time,
        d_ave_sim_sym_ltt$dd,
        col = "#ffc20a", 
        lwd = 2)
  lines(d_ave_sim_sym_ltt$time,
        d_ave_sim_sym_ltt$pd,
        col = "#1a85ff", 
        lwd = 2)
  lines(d_ave_sim_sym_ltt$time,
        d_ave_sim_sym_ltt$vi,
        col = "black", 
        lwd = 2)
  
  # add axes
  axis(
    side = 1,
    at = c(-190, -150, -100, -50, 0),
    col = "black",
    col.axis = "black",
    pos = round(min(d_ave_sim_sym_ltt$dd), 1) - 1.5,
    las = 1,
    cex.axis = 0.7,
    padj = -2,
    tck = -0.005
  )
  axis(
    side = 4,
    seq(round(min(d_ave_sim_sym_ltt$dd), 1), round(max(d_ave_sim_sym_ltt$dd), 1) + 6, by = 10),
    col = "black",
    col.axis = "black",
    pos = round(min(d_ave_sim_sym_ltt$dd), 1) + 2,
    las = 2,
    cex.axis = 0.7,
    hadj = 0.5,
    tck = -0.005
  )
  
  # add axes labels
  mtext("Time",
        side = 1,  
        line = 1,
        cex = 0.8)
  mtext("Lineages",
        side = 4, 
        line = 1,
        cex = 0.8)
  
  # add a legend
  legend("topleft",
         legend = c("Biphasic", "Direct development", "Paedomorphic", "Viviparous"),
         bty = "n",
         pch = 22,
         pt.bg = cols,
         pt.cex = 1.2,
         cex = 0.7)
}




