library(tidyverse)

setwd("C:/Users/eKrys/Downloads/")
lizards <- read.csv("lizards.csv")

lizards$time <- factor(lizards$time, levels = c("early", "midday", "late") )

#full interaction model
#  The effect of light
#  The effect of time
#  The effect of light and time together (light:time)
lm_int <- lm(grahami~light*time, data = lizards)

#This will keep the left side of the ~ (i.e., grahami) and 
# replace the right side with light + time, instead of light*time
#   The effect of light (alone)
#   The effect of time (alone)
lm_add <- update(lm_int, . ~ light + time)

#Taking lm_add, and removing time
#   The effect of light alone
lm_light <-update(lm_add, . ~ . -time)

#Taking lm_add, and removing light
#   The effect of time alone
lm_time <- update(lm_add, . ~ . -light)

lm_null <- update(lm_add, . ~ 1) #1 is symbolic for the intercept

#Look at the summary statistics for the additive model

summary(lm_add)
summary(lm_null)
#This should be measuring the effect of light or time on grahami
#Intercept is the baseline effect; our baseline is shady + early
#lightsunny is measuring the effect of sunny conditions holding time constant (baseline, i.e., early)
#timemidday is measuring the effect of midday time holding light constant (baseline, i.e., shady)
#p-values provide statistical clarity that the effect is not zero
#This linear model is:  
#grahami = 27.288 + -19.325(sunny) + 13.137(midday) + -9.500(late)

#Test time by comparing with light model (drop time)
#ANOVA for regression is not the same as ANOVA for groups; here the question is
#Does adding the predictor, time, improve the model enough to justify it's 
# extra complexity?

#We are comparing grahami ~ light + time, to grahami ~ light
#ANOVA will compare:
#       How much residual sum of squares (RSS) decreases
#       relative to how many extra parameters you added
#       relative to the remaining noise in the model
#Test two models against each other; 
# should test two nested models; where one has some parameters that the other
# does not, to check if one is better...
anova(lm_add, lm_light)
car::Anova(lm_add, lm_light)

#Main things to look at;
#    F stats are not that useful on their own
#    The point of this is the Df (2 parameters differencces between the models)
#    in Summary(lm_add) I get early vs. midday and early vs. late; but what about
#    midday vs. late;
#    When you have lots of levels you have lots of comparisons; looking at each 
#    pairwise is bad since it is multiple comparisons
#   The ANOVA gives one test for the overall effect of time
#Results give a p-value of 0.01698, which is < 0.05
# Res.DF  Time adds 2 dummy variables, so you lose 2 df
# RSS     Residual sum of squares is much lower with time added
#           by 1912.3, which means a better fit 
# Df      Model 1 has 2 more parameters than model 2 (midday and late dummy vars)
# Pr(>F)  p-value; here it is < 0.05 so statistically clear decrease in RSS

#Testing light (identical to test in summary())

anova(lm_add, lm_time)

##In this simple case, drop1() and car::Anova will do exactly the 
## same thing we did above...

drop1(lm_add, test = "F")
drop1(lm_int, test = "F")
#Drop any single terms in the model;
#Try dropping each term, everytime going back to the full model;
#Results indicate that both time, and light make significant contributions to 
# the model and should be kept. 

#Better ANOVA table

car::Anova(lm_add)
car::Anova(lm_int)
#There is an effect of light, time, and light:time
#It's doing type II tests;
# Tests the interaction, then the effect of light:time against the additive model
# 

#Best thing to do when doing type III tests:
#    contrasts = c("contr.sum", "contr.poly"))

#########
#Interaction Model

print(summary(lm_int))

#The intercept is the baseline; shady-early
#The sunny is the differences between baseline holding time constant... 
#

#The p-values indicate if the difference between the parameter effect and the 
# baseline level is 0. For timemidday there is not statistically clear evidence
# that lightsunnyhas an effect on grahami, compared to baseline
# In contrast, timemidday appears to have a statistically clear effect on grahami
# compared to baseline (shady/early)

#Interaction terms
#Does the effect of light depend on time?
#Does the effect of time depend on light?

#lightsunny:timemidday -- does the effect of sunny vs. baseline light change
# when you move from baseline time to midday?
#lightsunny:timemidday = -32.833, p = 0.09266
#    This means "the difference between sunny and baseline light is 32.8 units
#                 more negative at midday than at baseline time:                  

#lightsunny:timelate -- does the effect of sunny vs. baseline light change
# when you move from baseline time to late?


#The full model is: grahami_i = β0 + β1 * lightsunny_i
#                                  + β2 * timemidday_i
#                                  + β3 * timelate_i
#                                  + β4 * (lightsunny_i * timemidday_i)
#                                  + β5 * (lightsunny_i * timelate_i)
#                                  + ε_i

#The interaction variables are same dummy variables; if sunny and midday = 1,
#                                                    if sunny and baseline = 0

car::Anova(lm_int)
#All the variables are contributing to the model

######
#Contrasts
library(emmeans)

#Estimated marginal means -- (least-squares means)

#Asking for estimated marginal means (EMMs) of light, from a model containing
# main effect of light, 
# main effect of time, 
# interaction between light and time. 
#When an interaction exists, the effect of light depends on time. 
#emmeans decides "how do I summarize the effect of light when it varies across 
# time?". The default answer is to average over levels of time. 

em_light <- emmeans(lm_int, ~ light)
#emmeans is warning us that light is not constant over time; which we saw with the 
# lm_int, where baseline to sunny at midday was different to baseline to sunny at 
# the baseline level. 

print(em_contr <- contrast(em_light, "eff"))

#the "eff" contrast computes effect coding:
# It compares each level to the grand mean of the factor
# The effects sums to zero
#Since "light" has two levels, the effects are:
#    Shady effect = +10.3
#    Sunny Effect = -10.3
#These are not differences between shady and sunny.They are deviations 
# from the overall average of light, averaged across time.
#The p-values are testing: Is this level's mean significantly different from the 
# grand mean (averaged over time)?

#Why does this matter?
# The model has a significant interaction, the effect of light:
#   it is not constant over time
#   it is not well summarized by a single number
#   should be ideally examined at each time level

plot(emmeans(lm_int, ~ light))

plot(em_contr)


#Effects
#A main effect is the effect of a factor when you hold other factors constant
# e.g., the effect of light, or the effect of time, when other factors are 
# being held constant. 

#If an interaction is present (i.e., time:light) then main effects are conditional
# The main effects become effect of light when time = baseline, the effect of time
# when light = baseline. They are no longer overall effects of light or time, 
#  they are simple effects at the baseline level of other variables. 

#Treatment effects (ANOVA or Regression)
# Treatment effects are the differences between the mean of a particular group 
# and some reference value (baseline, grand mean, another group). 

#Examples of treatment effects:
#    shady vs. sunny
#    midday vs early
#    late vs. early
#     Shady-midday vs. Shady-early
#These are comparisons between groups. 

#Treatment effecsts are what you test when you run pairwise comparisons or 
# contrasts, emmeans comparisons, or ANOVA F-tests for factors. These answer
# questions like "Is midday sig. different from early?", "Is sunny sig. different from 
# shady?". "Is the efefct of light different at midday than early?".

#With our model (lm_int) there is a significant interaction between light and time
# This means main effects are relative to baseline conditions, and these are not 
# overall effects. 

#Treatment effects come from emmeans() and contrasts

#Give me the effect of light at each time level, using the fitted interaction model
plot(emmeans(lm_int, ~light | time))

#RESULTS

time = early:
light emmean   SE df lower.CL upper.CL
shady  23.50 5.38 17   12.159     34.8
sunny  11.75 5.38 17    0.409     23.1

#At time = early, shady is 11.75 units higher than sunny
#this matches the effect of lightsunny because early is the baseline
#At early time, sunny conditions reduces grahami moderately
#Based on confidence intervals this is probably statistically unclear

time = midday:
light emmean   SE df lower.CL upper.CL
shady  51.25 5.38 17   39.909     62.6
sunny   6.67 6.21 17   -6.429     19.8

#At time = midday, shady is 44.6 units higher than sunny
# This huge difference is expected from the interaction term
# (i.e., lightsunny:timemidday = -32.833)
#At midday, sunny light results in lower grahami numbers
#This effect is probably statistically clear as well

time = late:
light emmean   SE df lower.CL upper.CL
shady  10.75 5.38 17   -0.591     22.1
sunny   5.50 5.38 17   -5.841     16.8

#At time = late, shady is 5.25 units higher than sunny. Small difference is
# consistent with the other interaction term (i.e., lightsunny:timelate = +6.5). 
#At late time, an effect of light is statistically unclear and small


#Doing the stats
pairs(emmeans(lm_int, ~ light | time))

#Only sunny conditions at midday is statistically clear, so our conclusions are
# It is statistically clear that grahami observations are lower at midday in sunny
# conditions. 

em <- emmeans(lm_int, ~ light | time)
print(em)

library(ggplot2)

em_df <- as.data.frame(emmeans(lm_int, ~ light | time))

ggplot(em_df, aes(x = time, y = emmean, color = light, group = light)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15) +
  labs(
    x = "Time of Day",
    y = "Predicted grahami",
    color = "Light Condition",
    title = "Interaction Plot: Light × Time"
  ) +
  theme_minimal(base_size = 14)


install.packages("PlantGrowth")
library(PlantGrowth)

data <- PlantGrowth
data