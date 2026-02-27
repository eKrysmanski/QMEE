library(tidyverse)

#Reading in the data...
BQ_compile <- read.csv("data/SAR-assay-data/combined-flg22-data.csv")

#Setting genotype as a factor with an order vis. levels = 
BQ_compile$athal_genotype <- factor(BQ_compile$athal_genotype)
BQ_compile$treatment <- factor(BQ_compile$treatment, levels = c("MgCl2", 
                                                                "Pst AvrRpt2", 
                                                                "H2O", 
                                                                "Flg22"))
BQ_compile$population <- as.numeric(BQ_compile$population)
BQ_compile$induction <- factor(BQ_compile$induction, levels = c("Mock", "SAR"))

str(BQ_compile)

#Creating an ID for each based on combination of genotype and treatment
BQ_compile <- BQ_compile %>% 
  mutate(geno_treat = paste(athal_genotype, treatment, sep = "_"))

#Take a look at the data
library(ggplot2)

ggplot(data = BQ_compile, aes(x = treatment, y = population, color = treatment)) +
  geom_point()

#Generate the glm...
#This is basically count data, counts of colony forming units per leaf disk
#A poisson distribution should be appropriate for this kind of data; using the 
# quasipoisson since Dr. Bolker said there is no negatives to using it, and it 
# can help account for overdispersion (so might as well default to it)
library(performance)

BQ_glm <- glm(population ~ treatment, 
              data = BQ_compile, 
              family = quasipoisson(link = "log"))

BQ_glm2 <- update(BQ_glm, . ~ treatment + experiment)

anova(BQ_glm, BQ_glm2, test = "F")

#Try a negative binomial family to see what it is like
library(MASS)

BQ_glm_nb <- glm.nb(population ~ treatment, 
                    data = BQ_compile)

BQ_glm_nb2 <- update(BQ_glm_nb, . ~ treatment + experiment)

#Diagnostic plots
performance::check_model(BQ_glm2)
performance::check_model(BQ_glm_nb2)

plot(BQ_glm2)
plot(BQ_glm_nb2)

par(mfrow=c(1,2))
plot(BQ_glm2, which=3, main="QP: Residuals vs Fitted")
plot(BQ_glm_nb2, which=3, main="NB: Residuals vs Fitted")


summary(BQ_glm2)
summary(BQ_glm_nb2)


#Inferential Plots

dotwhisker::dwplot(list(QP_model = BQ_glm2, 
                        NB_model = BQ_glm_nb2))

dotwhisker::dwplot(BQ_glm_nb2)
dotwhisker::dwplot(BQ_glm2)

library(emmeans)
library(ggrepel)

plot(emmeans(BQ_glm_nb2, ~ treatment + experiment, comparison = TRUE))

emms <- emmeans(BQ_glm_nb2, ~ treatment, level = 0.834)

plot_data <- as.data.frame(emms)

ggplot(plot_data, aes(x = treatment, 
                      y = emmean)) +
  geom_point(size = 2.5, 
             position = position_dodge(width = 0.25)
  ) +
  geom_line(linewidth = 1, 
            position = position_dodge(width = 0.25),
            alpha = 0.75
  )+
  geom_errorbar(aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), 
                position = position_dodge(width = 0.25), 
                linewidth = 0.75,
                width = 0.25,
                alpha = 0.75
  ) +
  ylab(label = "Predicted Population [log(cfu/ld)]") +
  xlab(label = "Treatment") +
  labs(title = "Predicted Means [log10(cfu/ld)]") +
  geom_text_repel(aes(label = c("MgCl2", "Pst AvrRpt2", "H20", "Flg22")), 
                  nudge_x = 0.3,
  ) +
  theme_bw() +
  theme(
    legend.position = "none", 
    
  )

dotwhisker::dwplot(BQ_glm_nb2, 
                   dot_args = list(size = 2.5, 
                                   shape = 18, 
                                   color = "darkred"), 
                   whisker_args = list(size = 4, 
                                       alpha = 0.5)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(title = "Coefficient Plot") +
  theme(
    legend.position = "none"
  )
