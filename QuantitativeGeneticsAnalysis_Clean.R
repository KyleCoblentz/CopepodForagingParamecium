################################################################################
### Examine heritability and trait correlations
################################################################################

### load packages

library(tidyverse); library(ggcorrplot); library(cowplot)

### load phenotype data

phenotype <- read.csv('Paramecium_Phenotype_Data.csv')

### look at correlation matrix 

# drop variables we don't want first

str(phenotype)

phenotype <- phenotype %>% select(-c(X, mean_grey, sd_grey,
                                     sd_area, sd_perimeter,
                                     sd_major, sd_minor,
                                     sd_ar, duration, N_frames, id))


phenotype_median <- phenotype %>% group_by(file) %>% summarise_all(.funs = function(x) median(x))

corr_median <- cor(phenotype_median[,-1])

median_corr_matrix_plot <- ggcorrplot(corr = corr_median)

# save_plot(filename = 'corr_matrix.png',plot = mean_corr_matrix_plot,
#          bg = 'white', nrow = 1.35, ncol = 1.35)

### keep sd_gross speed, gross_speed, net_disp, 
### sd_turning, mean_turning, and all of the size variables

phenotype <- phenotype %>% select(c(file, mean_major,
                                    mean_minor, mean_ar, mean_turning,
                                    sd_turning, gross_speed, net_disp,
                                    sd_gross_speed))

# create a column that is genotype without all of the front information

phenotype <- phenotype %>% mutate(genotype = substring(file, first = 11))

### ANOVA's to look at heritability of traits

### will do on the phenotype data including all individuals

### length of paramecium

fit_major <- lm(mean_major ~ genotype, data = phenotype)

anova_fit_major <- anova(fit_major)

anova_fit_major$`Sum Sq`[1]/sum(anova_fit_major$`Sum Sq`)

# h2 = 0.62

### width

fit_minor <- lm(mean_minor ~ genotype, data = phenotype)

anova_fit_minor <- anova(fit_minor)

anova_fit_minor$`Sum Sq`[1]/sum(anova_fit_minor$`Sum Sq`)

# h2 = 0.59

### aspect ratio

fit_ar <- lm(mean_ar ~ genotype, data = phenotype)

anova_fit_ar <- anova(fit_ar)

anova_fit_ar$`Sum Sq`[1]/sum(anova_fit_ar$`Sum Sq`)

# h2 = 0.56

### mean turning

fit_turning <- lm(mean_turning ~ genotype, data = phenotype, )

anova_fit_turning <- anova(fit_turning)

anova_fit_turning$`Sum Sq`[1]/sum(anova_fit_turning$`Sum Sq`)

# h2 = 0.14

### sd_turning

fit_sd_turning <- lm(sd_turning ~ genotype, data = phenotype, )

anova_fit_sd_turning <- anova(fit_sd_turning)

anova_fit_sd_turning$`Sum Sq`[1]/sum(anova_fit_sd_turning$`Sum Sq`)

# h2 = 0.38

### gross_speed

fit_gross_speed <- lm(gross_speed ~ genotype, data = phenotype)

anova_fit_gross_speed <- anova(fit_gross_speed)

anova_fit_gross_speed$`Sum Sq`[1]/sum(anova_fit_gross_speed$`Sum Sq`)

# h2 = 0.72

### net displacement

fit_net_disp <- lm(net_disp ~ genotype, data = phenotype)

anova_fit_net_disp <- anova(fit_net_disp)

anova_fit_net_disp$`Sum Sq`[1]/sum(anova_fit_net_disp$`Sum Sq`)

# h2 = 0.14


### sd gross speed

fit_sd_gross_speed <- lm(sd_gross_speed ~ genotype, data = phenotype)

anova_fit_sd_gross_speed <- anova(fit_sd_gross_speed)

anova_fit_sd_gross_speed$`Sum Sq`[1]/sum(anova_fit_sd_gross_speed$`Sum Sq`)

### h2 = 0.61

