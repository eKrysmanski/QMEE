#Assignment 6 - Linear Models
library(tidyverse)
library(performance)


data <- read.csv("fmo1-ald1-pooled.csv", header = TRUE)

data$genotype <- factor(data$genotype, levels = c("col-0", "fmo1-1", "ald1-T2"))

lm_int <- lm(population ~ genotype*treatment, data = data)

#Check diagnostics
performance::check_model(lm_int, check = c("homogeneity"))

#Not great, not bad...homogeneity of variance may not be great but the 
# data I have is very right skewed since there are two SAR defective phenotypes 
# with the 3 mutants. It's also bacterial level data so it would very likely 
# benefit from a log transformation

#Check the base R versions since im more used to these
plot(lm_int, which = 1) #Funnel shape; looks like increasing variance (heteroscedasticity)
plot(lm_int, which = 2) #Normal-ish; looks like a right skew? 
plot(lm_int, which = 3) #Does not look good, increasing trend indicating 
                        # increased variance with higher mean... 

#Bacterial data often get's log transformed since it is exponential...

log_data <- data %>% 
  mutate(population = log10(population))

lm_int <- lm(population ~ treatment*genotype, data = log_data)

#Check diagnostics
performance::check_model(lm_int, check = c("homogeneity")) #looks ok I think...

#Checking base R diag plots:
plot(lm_int, which = 1) #much less funnel shaped; good enough I think (homoscedastic)
plot(lm_int, which = 2) #Normal-enough, minor skew but normality matters the least 
plot(lm_int, which = 3) #very slight increasing trend; but not bad enough to worry I think

plot(lm_int) #Nothing really sitcks out with constant leverage... 

#Check model summary and explain the results:
summary(lm_int)

#Note: the p value is for  "Is this coefficient clearly different from zero?" 

#Looks like treatmentSAR results in reduction in bacterial level (makes sense), this 
# reduction is statistically clear based on the p-value. 
#Looks like it is not statistically clear if the coefficients for genotypeald1-T2 
# and genotypefmo1-1 are different from zero. Infering from the mock treatments 
# that the infection represents local defense, this is not what I predicted. Both
# of these genotypes have a small but statistically clear defect in local defense. 
# Based on my data, this is not observed. However, this assay is a SAR assay, and
# not designed for assessing defects in local defense, so it must be taken with a 
# grain of salt. There does appear to be a small positive effect, but it is not 
# statistically clear. 
#Looks like treatmentSAR:genotypefmo1-1 and ald1-1 are counteracting the SAR reduction
# This makes sense since these genotypes are SAR defective, so the reduction in 
# the response variable from the treatmentSAR coefficient must be countered when
# the genotype is either fmo1-1 or ald1-T2. 

anova(lm_int)
#Looks like all the variables make significant contributions to the model;


#Inferential plot
library(dotwhisker)

dwplot(lm_int, 
       dot_args = list(size = 2.5, 
                       shape = 18, 
                       color = "darkred"), 
       whisker_args = list(size = 4, 
                           alpha = 0.5)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  xlab(label = "Coefficient Values Log10(cfu/ld)") +
  labs(title = "Coefficient Plot") +
  theme(
    legend.position = "none"
  )

#Prediction Plot
library(emmeans)

plot(emmeans(lm_int, ~genotype*treatment), comparison = TRUE)

#This plot makes sense; but I think there is a better way to visualize this info. 
#Since I'm most interested in changes in the response variable between mock/SAR
# for each genotype, try and plot each genotype in a lineplot with mock and SAR
# as the independent variables... 

library(ggrepel)

#Better? Prediction Plot
emmeans <- emmeans(lm_int, ~ genotype * treatment)
plot_data <- as.data.frame(emmeans)

ggplot(plot_data, aes(x = treatment, 
                      y = emmean, 
                      color = genotype, 
                      group = genotype)) +
  geom_point(size = 2.5, 
             position = position_dodge(width = 0.25)
             ) +
  geom_line(linewidth = 1, 
            position = position_dodge(width = 0.25),
            alpha = 0.75
            )+
  geom_errorbar(aes(x = treatment, ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(width = 0.25), 
                linewidth = 0.75,
                width = 0.25,
                alpha = 0.75
                ) +
  ylab(label = "Predicted Population [log(cfu/ld)]") +
  xlab(label = "Treatment") +
  labs(title = "Predicted Means [log10(cfu/ld)]") +
  geom_text_repel(aes(label = c("", "", "", "Col-0", "fmo1-1", "ald1-T2")), 
                   nudge_x = 0.3,
                   ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )


#Looks like genotypefmo1-1 and genotype ald1-T2 have very little effect on the 
#response variable; they are quite close to 0 with CIs crossing 0. This suggests 
#a small but statistically unclear effect on genotype. This result is slightly 
#surprising, as there are reports of a weak defect for local defense responses for
#these mutants. However, this assay is not testing local defense defects, but 
#SAR responses, so I cannot make any conclusions biologically about this. 
#Looks like treatmentSAR has a strong negative effect on the response variable; 
#this makes sense since the baseline is Col-0 mock, and holding so SAR treatment
#should have a strong negative effect on the bacterial level in distant leaves. 
#The interaction variables, treatmentSAR:genotypefmo1-1/ald1-T2 both have strong
#positive effects. This appears to be to counter the effect of SAR; This indicates
#that SAR has no effect on the bacterial level for these plants, which aggrees with 
#observations in the literature. 

#Looks like there is a right skew, and

#Quick peek at the data  
data_sum <- data %>% 
  mutate(ID = paste(genotype, treatment)) %>% 
  group_by(ID) %>% 
  summarize(mean = mean(population)) 

data_sum$ID <- factor(data_sum$ID, levels = c("col-0 mock", "col-0 SAR", 
                    "fmo1-1 mock", "fmo1-1 SAR", 
                    "ald1-T2 mock", "ald1-T2 SAR"))

ggplot(data = data_sum, aes(x = ID, y = mean)) +
  geom_col() +
  scale_y_log10()

  