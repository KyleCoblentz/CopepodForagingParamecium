################################################################################
### Copepod Foraging Analysis
################################################################################

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(brms); library(cmdstanr); library(tidyr); library(ggbiplot)

### load data

feed <- read.csv('Foraging trial data.csv')
  
phenotype <- read.csv('Paramecium_Phenotype_Data.csv')

copepod <- read.csv('CopepodForagingLengths_Clean.csv')

### first focus on copepod foraging data

# make prey type a factor

feed$Prey.type <- as.factor(feed$Prey.type)

# get rid of observation with dead copepod

feed <- feed %>% filter(Dead == 0)

# create a proportion eaten variable

feed <- feed %>% mutate(Prop.Eaten = 1 - Prey.remaining/Prey.offered)

# change Prey.type to genotype

colnames(feed)[4] <- 'genotype'

################################################################################
### match phenotype data to foraging data
################################################################################

### perform pca on the medians of phenotypes across all of the genotypes

# drop columns of phenotype we don't want

phenotype <- phenotype %>% select(-c(X, mean_grey, sd_grey, sd_area, sd_perimeter, 
                                     sd_major, sd_minor, sd_ar, duration, N_frames,
                                     id))

### select the remaining phenotypes that we do want

phenotype <- phenotype %>% select(c(file, mean_major, mean_minor, mean_ar, mean_turning,
                                    sd_turning, gross_speed, sd_gross_speed, 
                                    net_disp))

# modify phenotype data to be medians

phenotype <- phenotype %>% group_by(file) %>% summarise_all(list(median))

# create a column that is genotype without all of the front text

phenotype <- phenotype %>% mutate(genotype = substring(file, first = 11))

# run a pca on the phenotype data

pca_pheno_data <- phenotype %>% select(-c(file, genotype))

pheno_pca <- prcomp(pca_pheno_data, center = TRUE, scale = TRUE)

summary(pheno_pca)


### create a pca plot

# make a dataframe for pca 

# we will take the negative of pca's 1 and 2 so that they can be more easily 
# interpreted as increasing speed and increasing size

pca_plot_data <- data.frame(pca1 = -pheno_pca$x[,1], pca2 = -pheno_pca$x[,2], 
                            in_foraging = ifelse(phenotype$genotype %in% feed$genotype, 1, 0), 
                            genotype = phenotype$genotype)

# want to plot the points as the genotype names and as different colors depending
# on whether they were included in the analysis or not

pca_base <- ggplot(data = pca_plot_data, aes(x = pca1, y = pca2, color = as.factor(in_foraging), shape = as.factor(in_foraging))) + 
  geom_point(size = 2) + scale_color_manual(values = c('1' = 'magenta', '0' = 'black'),
                                    labels = c('no', 'yes'), name = 'In foraging trial?') + 
  scale_shape_manual(values = c('1' = 15, '0' = 16), labels = c('no', 'yes'), name = 'In foraging trial?') +
  theme_cowplot() + 
  xlab('PC 1 \n(37.2% explained variance)') + ylab('PC 2 \n(23% explained variance)')
  
### want to add the axis loadings as vectors

# set up what we need by manipulating code from the package ggbiplot

# from ggbiplot

nobs.factor <- sqrt(nrow(pheno_pca$x) - 1)

d <- pheno_pca$sdev

u <- sweep(pheno_pca$x, 2, 1/(d * nobs.factor), FUN = '*')

v <- pheno_pca$rotation

choices <- pmin(1:2, ncol(u))

df.u <- as.data.frame(sweep(u[,choices], 2, d[choices], FUN = '*'))

df.v <- as.data.frame(v[,choices])

names(df.v) <- c('xvar', 'yvar')

names(df.u) <- c('xvar', 'yvar')

df.u <- df.u * nobs.factor

r <- sqrt(qchisq(0.69, df = 2)) * prod(colMeans(df.u^2))^(1/4)

v.scale <- rowSums(v^2)

df.v <- r*df.v/sqrt(max(v.scale))

df.v$varname <- rownames(v)

df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
df.v$hjust <- with(df.v, (0.25 * sign(xvar))/2)

### the following code will adjust how variables names are plotted onto 
### figure

df.v$hjust[c(2,4)] <- 1.1

### adjust the locations of the variable names for the pca loadings

df.v$xvar_mod <- 0.8*-df.v$xvar + c(0,0,0,0,0,0,0,0.35)

df.v$yvar_mod <- 0.8*-df.v$yvar + c(0, 0, 0.6, 0, 0, 0, 0.4, -0.45)

df.v$angle_mod <- c(df.v$angle[1:2], 0 , df.v$angle[4:5], 0, 0, 0)

### put all of this together to create the plot with axis loading vectors

pca_axes_plot <- ggplot(data = pca_plot_data, aes(x = pca1, y = pca2, color = as.factor(in_foraging))) + 
  geom_point(aes(shape = as.factor(in_foraging)), alpha = 0.1, size = 2) + scale_color_manual(values = c('1' = 'magenta', '0' = 'black'),
                                    labels = c('no', 'yes'), name = 'In foraging trial?') + 
  scale_shape_manual(values = c('1' = 15, '0' = 16), labels = c('no', 'yes'), name = 'In foraging trial?') +
  theme_cowplot() + 
 geom_segment(data = df.v, aes(x = 0, y = 0, xend = 0.8*-xvar, yend = 0.8*-yvar),
                                         arrow = arrow(length = unit(1/2, 'picas')),
                                         color = muted('red')) + 
  geom_text(data = df.v, aes(label = varname, x = xvar_mod,
                             y = yvar_mod, angle = angle_mod, hjust = hjust), inherit.aes = FALSE) + 
  xlab('PC 1 \n(37.2% explained variance)') + ylab('PC 2 \n(23% explained variance)') + xlim(c(-6,5))

### put the two pca plots together

pca_plots <- plot_grid(pca_base, pca_axes_plot, ncol = 1, labels = 'AUTO')

 save_plot(filename = 'pca_plots.png', plot = pca_plots, nrow = 2, ncol = 1,
          bg = 'white', base_height = 3.5, base_width = 6.3)


################################################################################
### put phenotype and feeding trial data together for analysis
################################################################################

### extract the pca values for the first two axes 

phenotype <- phenotype %>% mutate(pca1 = -pheno_pca$x[,1], pca2 = -pheno_pca$x[,2]) 

###add the pca values to the feed data frame

phenotype_pca_merge <- phenotype %>% select(genotype, pca1, pca2)

feed <- left_join(feed, select(phenotype, c(pca1, pca2, genotype)), by = 'genotype')

### create an eaten column

feed$Prey.eaten <- feed$Prey.offered - feed$Prey.remaining

### add the copepod lengths to the analysis

feed <- left_join(feed, copepod, by = 'Copepod.ID')

### now we can fit a model 

### specify model priors

fit_prior <- c(prior(cauchy(0, 2), class = 'b'),
               prior(cauchy(0, 2), class = 'Intercept'),
               prior(cauchy(0, 2), class = 'sd'))

# fit model

fit <- brm(formula = Prey.eaten | trials(40) ~ 0 + 1 + pca1*pca2 + Length + (1|genotype) + (1|Copepod.ID),
           data = feed, family = beta_binomial(link = 'logit', link_phi = 'log'), prior = fit_prior, cores = 4, iter = 4000, warmup = 2000,
           backend = 'cmdstanr')

prior_summary(fit)

summary(fit, prob = 0.90, robust = TRUE)

fixef(fit, probs = c(0.05, 0.95), robust = TRUE)

### construct plots 

# predictions for pca1

# get predictions with pca2 set to different levels

# what values for pca2

quantile(feed$pca2, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

# setup new data

pca1_newdat <- data.frame(pca1 = rep(seq(-6.5, 3.4, by = 0.05), times = 3),
                          pca2 = c(rep(-1.9, times = 199), rep(-0.15, times = 199), rep(1.05, times = 199)),
                          pca2_quantile = c(rep('0.1', times = 199), rep('0.5', times = 199), rep('0.9', times = 199)),
                          Length = rep(mean(feed$Length, na.rm = TRUE), times = 597))

pred_pca1 <- fitted(fit, newdata = pca1_newdat, re_formula = NA, probs = c(0.05, 0.95))

pred_pca1 <- cbind(pred_pca1, pca1_newdat)

feed_plot <- feed %>% drop_na(pca1, Length)

feed_plot$pca2_quantile <- ifelse(feed_plot$pca2 < quantile(feed_plot$pca2, probs = 0.15), '0.1',
                                  ifelse(feed_plot$pca2 > quantile(feed_plot$pca2, probs = 0.85), '0.9', '0.5'))  

PCA1_plot <- ggplot(data = pred_pca1, aes(x = pca1, y = Estimate/40, color = pca2_quantile, shape = pca2_quantile)) + 
  geom_line() + geom_ribbon(aes(ymin = Q5/40, ymax = Q95/40, fill = pca2_quantile), alpha = 0.25) + 
  theme_cowplot() + geom_point(data = feed_plot, aes(x = pca1, y = Prey.eaten/40), size = 2)+ scale_color_brewer(palette = 'Dark2', labels = c('0.1 (Small)', '0.5', '0.9 (Large)')) + 
  scale_fill_brewer(palette = 'Dark2', labels = c('0.1 (Small)', '0.5', '0.9 (Large)')) + 
  scale_shape_manual(values = c('0.1' = 16, '0.5' = 17, '0.9' = 18), labels = c('0.1 (Small)', '0.5', '0.9 (Large)')) + guides(fill = guide_legend(title = 'PCA2 (Size) Quantile'), color = guide_legend(title = 'PCA2 (Size) Quantile'),
                                                shape = guide_legend(title = 'PCA2 (Size) Quantile')) +
  xlab('Principal Component Axis (PCA) 1 (Speed)') + ylab('Proportion Paramecium Eaten')

# predictions for pca2

# get predictions with pca1 set to different levels

# what values for pca1

quantile(feed$pca1, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

# setup new data

pca2_newdat <- data.frame(pca2 = rep(seq(-4, 4.6, by = 0.05), times = 3),
                          pca1 = c(rep(-1.85, times = 173), rep(0.15, times = 173), rep(1.71, times = 173)),
                          pca1_quantile = c(rep('0.1', times = 173), rep('0.5', times = 173), rep('0.9', times = 173)),
                          Length = rep(mean(feed$Length, na.rm = TRUE), times = 173*3))

pred_pca2 <- fitted(fit, newdata = pca2_newdat, re_formula = NA, probs = c(0.05, 0.95))

pred_pca2 <- cbind(pred_pca2, pca2_newdat)

feed_plot$pca1_quantile <- ifelse(feed_plot$pca1 < quantile(feed_plot$pca1, probs = 0.15), '0.1',
                                  ifelse(feed_plot$pca1 > quantile(feed_plot$pca2, probs = 0.85), '0.9', 0.5))  


PCA2_plot <- ggplot(data = pred_pca2, aes(x = pca2, y = Estimate/40, color = pca1_quantile, shape = pca1_quantile)) + 
  geom_line() + geom_ribbon(aes(ymin = Q5/40, ymax = Q95/40, fill = pca1_quantile), alpha = 0.25) + 
  theme_cowplot() + geom_point(data = feed_plot, aes(x = pca2, y = Prey.eaten/40), size = 2) + scale_color_brewer(palette = 'Dark2', labels = c('0.1 (Slow)', '0.5', '0.9 (Fast)')) + 
  scale_fill_brewer(palette = 'Dark2', labels = c('0.1 (Slow)', '0.5', '0.9 (Fast)')) + 
  scale_shape_manual(values = c('0.1' = 16, '0.5' = 17, '0.9' = 18), labels = c('0.1 (Slow)', '0.5', '0.9 (Fast)')) + 
  guides(fill = guide_legend(title = 'PCA1 (Speed) Quantile'), color = guide_legend(title = 'PCA1 (Speed) Quantile'),
                                                                                               shape = guide_legend(title = 'PCA1 (Speed) Quantile')) + 
  xlab('Principal Component Axis (PCA) 2 (Size)') + ylab('Proportion Paramecium Eaten')

### posterior predictive check -- does predicted data look like the observed data?

post_check <- pp_check(fit, type = 'dens_overlay', ndraws = 100)

post_check <- post_check + xlab('Number of Prey Eaten') + ylab('Density') + 
  scale_color_manual(name = '',labels = c('Observed', 'Posterior Prediction\nof Data'), values = c('darkblue', 'lightblue'))

# save_plot(filename = 'post_pred_check.png', plot = post_check)

### put plots together

# response plot

response_plot <- plot_grid(PCA1_plot, PCA2_plot, nrow = 1, ncol = 2, labels = 'AUTO')

save_plot('response_plot.png', response_plot, nrow = 1, ncol = 2, bg = 'white')

### calculations for table

samples <- as_draws_df(fit)

samples$chain <- rep(c('1', '2', '3', '4'), each = 2000)

samples$iteration <- rep(1:2000, times = 4)

sum(samples$b_Intercept < 0)/length(samples$b_Intercept)

sum(samples$b_pca1 > 0)/length(samples$b_pca1)

sum(samples$b_pca2 > 0)/length(samples$b_pca2)

sum(samples$`b_pca1:pca2` > 0)/length(samples$`b_pca1:pca2`)

sum(samples$b_Length > 0)/length(samples$b_Length)

### trace plots for all of the parameters 

trace_intercept <- ggplot(data = samples, aes(x = iteration, y = b_Intercept, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Intercept')

trace_slope_PCA1 <- ggplot(data = samples, aes(x = iteration, y = b_pca1, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Slope PCA Axis 1')

trace_slope_PCA2 <- ggplot(data = samples, aes(x = iteration, y = b_pca2, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Slope PCA Axis 2')

trace_slope_Interaction <- ggplot(data = samples, aes(x = iteration, y = `b_pca1:pca2`, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Interaction PCA1:PCA2')

trace_copepod_sd <- ggplot(data = samples, aes(x = iteration, y = sd_Copepod.ID__Intercept, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Standard Deviation \nCopepod Random Intercept')

trace_line_sd <- ggplot(data = samples, aes(x = iteration, y = sd_genotype__Intercept, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Standard Deviation \nOutcrossed Line Random Intercept')

trace_phi <- ggplot(data = samples, aes(x = iteration, y = phi, color = chain)) + geom_line() + 
  theme_cowplot() + xlab('Iteration') + ylab('Beta-binomial phi')

### put trace plots together

trace_plot <- plot_grid(trace_intercept, trace_slope_PCA1, trace_slope_PCA2, trace_slope_Interaction,
                        trace_copepod_sd, trace_line_sd, trace_phi, ncol = 2, nrow = 4)

# save_plot(filename = 'trace_plot.png', plot = trace_plot,
#          ncol = 2, nrow = 3, bg = 'white')

### plot of residuals versus pca2 for supplementary material.

# data for residual plot

feed_residual <- feed %>% filter(!is.na(Length)) %>% filter(!is.na(pca1))

residuals <- residuals(fit)

feed_residual <- feed_residual %>% mutate(Residuals = residuals[,1])

pca2_resid_plot <- ggplot(feed_residual, aes(x = pca2, y = Residuals)) + geom_point() + 
  xlab('PCA2 (Size)') + ylab('Average Residual') + theme_cowplot() + 
  geom_hline(yintercept = 0)

# save_plot(filename = 'residual_nonlinear_plot.png', plot = pca2_resid_plot)

### make a heatmap of the 'fitness' surface

max(feed$pca1, na.rm = TRUE)

min(feed$pca1, na.rm = TRUE)

# pca1 ranges from -5.7 to 3.4 

max(feed$pca2, na.rm = TRUE)

min(feed$pca2, na.rm = TRUE)

# pca2 ranges from -4.3 to 4.24 

# want to make a grid. For each value of pca1, replicate that value 
# 100 times. Then for pca 2 we will get each value 100 times to match 
# pca 1.

new_data_3D <- data.frame(pca1 = rep(seq(from = -5.7, to = 3.4, length.out = 100), each = 100),
                          pca2 = rep(seq(from = -4.3, to = 4.3, length.out = 100), times = 100),
                          Length = rep(mean(feed$Length, na.rm = TRUE), times = 100*100))

# now we need to get predictions from our model

predictions_3D <- fitted(fit, newdata = new_data_3D, re_formula = NA, probs = c(0.05, 0.95))

predictions_3D <- cbind(new_data_3D, predictions_3D)

# drop areas where there aren't data

predictions_3D$pca1 <- ifelse(predictions_3D$pca1 < -2 & predictions_3D$pca2 > 7.3 + predictions_3D$pca1*1.5, NA, predictions_3D$pca1)

predictions_3D$pca1 <- ifelse(predictions_3D$pca1 < -1.5 & predictions_3D$pca2 < -4.945 - predictions_3D$pca1*0.429, NA, predictions_3D$pca1)

predictions_3D$pca1 <- ifelse(predictions_3D$pca1 > 0 & predictions_3D$pca2 > 4.3 - predictions_3D$pca1*0.897, NA, predictions_3D$pca1)

predictions_3D$pca1 <- ifelse(predictions_3D$pca1 > 1 & predictions_3D$pca2 < -5.05 + predictions_3D$pca1*0.75, NA, predictions_3D$pca1)


fitness_surface <- ggplot(filter(predictions_3D, !is.na(pca1)), aes(x = pca1, y = pca2, fill = (40-Estimate)/40)) + geom_tile() + geom_contour(aes(z = (40-Estimate)/40)) + 
  scale_fill_viridis_c(name = 'Expected Fitness \n(Proportion Surviving)') + geom_point(data = feed, aes(x = pca1, y = pca2), inherit.aes = FALSE) + theme_cowplot() + 
  xlab('PCA1 (Speed)') + ylab('PCA2 (Size)') + theme(legend.title.align = 0.5)

save_plot(filename = 'Fitness_Surface.png', plot = fitness_surface, bg = 'white')


### plot of copepod size versus proportion of paramecium eaten for supplemental material

copepod_size_pred <- ggplot(data = feed, aes(x = Length, y = Prop.Eaten)) + geom_point() + 
  theme_cowplot() + xlab('Copepod Length') + ylab('Proportion Paramecium\nEaten')

# save_plot(filename = 'copepod_size_pred.png', plot = copepod_size_pred)


