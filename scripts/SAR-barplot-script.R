#### Barplot for SAR Assays -- Full Labs + Stats
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(showtext)
library(ggbeeswarm)

#Load in google font--------------------------------------------------------
font_add_google(name = "Pt Serif", family = "pt_serif") 
showtext_auto()

#Question: I really like Pt Serif as a font for figures, but I have read serif 
# fonts aren't great for figures - Is this just opinion or is there a real reason?

#Read in data
data <- read.csv("data/SAR-assay-data/SAR-flg22-2-2025-ANOVA.csv")

#Specify columns as factors with levels  
data$treatment <- factor(data$treatment, 
                         levels = c("MgCl2", "Pst AvrRpt2", "H2O", "Flg22"))

data$induction <- factor(data$induction, 
                         levels = c("Mock", "SAR"))

data$athal_genotype <- factor(data$athal_genotype, 
                              levels = c("Col-0"))

data$group <- factor(data$group, 
                     levels = c("1", "2"))

#Create a summarized data.frame 

data_sum <- data %>% 
  group_by(athal_genotype, treatment, group, induction) %>% 
  summarize(mean_pop = mean(population), 
            sd_pop = sd(population)) %>% 
  ungroup() %>% 
  as.data.frame()


#Inputs:________________________________________________________________________

#y-axis inputs--------------------------------
y_lims <- c(1e4, 5e7)
y_breaks <- 10^(4:7)

#Labels----------------------------------------
bar_labs <- c(
  expression(Mock),                  #Note: to stack use "expression(atop("text", MgCl[2]))
  expression(SAR),                   #Note: Use expression(italic("words")) for italics
  expression(Mock),                  #Note: Use expression(MgCl[2] to subscript
  expression(SAR)                    #Note: Use expression(words^superscript) for superscript
)

#group labels-----------------------------------
group_labs <- c(rep(expression(italic("Pst avrRpt2")), 2),   #This is a bit broken, for some reason needs a vector length same as data_sum 
                rep("Flg22 (1mM)", 2))                       #Works if you just specify each label identically twice since they stack

#Stats------------------------------------------
stats_labs <-c("a", "c", "a", "b") #One-way ANOVA + Tukey's HSD

#Adding labs to the data.frame (More robust than pulling labs from expression vectors)
data_sum$bar_labs <- bar_labs
data_sum$group_labs <- group_labs
data_sum$stats_labs <- stats_labs

################################################################################
############################ Plot Script #######################################

plot <- ggplot(data_sum, aes(x = group, y = mean_pop, fill = induction)) +
  geom_errorbar(                                               #Adding bar to represent mean
    position = position_dodge(0.8),                              #dodge
    aes(ymax = mean_pop,                                         #Set max to mean
        ymin = mean_pop),                                        #Set min to mean
    width = 0.3                                                  #Width of whisker
  ) +
  geom_errorbar(                                               #Adding standard deviation bars
    position = position_dodge(0.8),                              #Need to dodge
    aes(ymax = mean_pop + sd_pop,                                #Positive whiskers
        ymin = mean_pop - sd_pop),                               #Negative whiskers
    width = 0.2                                                  #width of the hats
    ) + 
  geom_beeswarm(                                               #Individual points jittered if overlapping
    data = data,                                                 #using raw dataframe
    aes(x = group, y = population, fill = induction),            #Same aesthetics as rest
    colour = "black",                                            #Outline colour
    shape = 21,                                                  #Circle with outline
    alpha = 0.7,                                                 #Transparency in case they overlap
    size = 4,                                                    #Point size
    stroke = 1,                                                  #border width
    dodge.width = 0.8                                            #geom_beeswarm does not accept position_dodge()
    ) +
  scale_y_log10(                                               #Adding a logarithmic scale
    limits = y_lims,                                             #Specify axis limits (set above)
    breaks = y_breaks,                                           #Specify label breaks (set above)
    labels = trans_format("log10", math_format(10^.x)),          #Adds superscript for labels
    oob = scales::oob_squish,                                    #Prevents out-of-bounds data from being removed
    expand =  expansion(0)                                       #Prevents cushion between bottom of bars and axis
    ) +
  geom_text(             #BAR LABELS                           #Adding bar-lables
    aes(y = 1,                                                   #Y-coordinate for label (Set to 1 so I can vjust below axis)
        label = bar_labs),                                       #The labels (set above)
    family = "pt_serif",                                         #Font
    size = 6.5,                                                  #Font size 
    position = position_dodge(0.8),                              #Need to dodge
    vjust = 1.5                                                  #vertical adjust location of text
  ) +
  geom_text(             #GROUP LABELS                         #Group labels
    aes(y = 1,                                                   #No position_dodge()
        label = group_labs),                                     #Specify in inputs   
    family = "pt_serif",                                         #Font
    size = 6.5,                                                  #Font size
    vjust = 3.25                                                 #Adjust below x-axis
  ) +
  geom_text(             #Label Bar                             #I acknowledge this is probably very stupid, 
    aes(y = 1,                                                   #but it works and is very easy to adjust as needed.  
        label = "_______________",                               #I have tried geom_hline(), annotate(), geom_segment()
        fontface = 2),
    family = "pt_serif", 
    size = 7, 
    vjust = 1.95
  ) + 
  geom_text(            #STAT LABELS                  #Adding in significance letters from Anova + Tukey's HSD    
    aes(y = mean_pop + sd_pop * 1.5,                    #Set y-coordinate at a fixed distance above each positive whisker
        label = stats_labs),                            #Labels are specified above
    family = "pt_serif",                                #Font
    size = 7,                                           #Font size
    position = position_dodge(0.8),                     #Dodge with the bars
    vjust = -1                                          #Adjust slightly higher
  ) +
  coord_cartesian(                  #This is required to make logticks work; I believe it basically zooms in
    clip = "off",                      #Allows for logticks to display outside
    expand = TRUE                      #Allows for space between the y-axis and bars
  ) +
  annotation_logticks(              #Log tick marks
    sides = "l",                       #On the left (l) side of plot
    outside = TRUE,                    #Ticks on outside - requires clip = "off" above
    linewidth = 0.9
  ) +
  scale_fill_manual(                                          #Create discrete scale for colour fill
    values = c("Mock" = "#35B749","SAR" = "#31688E")          #Supervisor really likes colours; she likes very bright and I prefer more mute...
    ) + 
  labs(
    y = expression("Bacterial density in distant leaves ( " ~ cfu ~ ld^-1*")")
  ) +
  theme_minimal(base_size = 14) +                  #Basic theme default
  theme(                                           #Adjusting theme features
    text = element_text(family = "pt_serif"),              #Font for labels
    panel.grid = element_blank(),                          #Remove grid lines
    axis.line = element_line(color = "black",              #Setting solid line for axes
                             linewidth = 0.9),             #Setting line width for axes
    axis.text.y = element_text(vjust = 0.2,                #Adjusting position and colour of y-labs
                               hjust = 2,                  #Horizontal adjustments
                               color = "black",            #Font colour
                               size = 20),                 #Font size
    axis.title.y = element_text(margin = margin(r = 15),   #Increase right margin
                                size = 20),                #Font size
    axis.title.x = element_blank(),                        #Remove the x-axis title,
    axis.text.x = element_blank(),                         #Remove x axis text
    plot.margin = margin(t = 20, r = 30, b = 50, l = 40),  #Set plot margins
    legend.position = "none"                               #Remove legend
  )

plot

#Curse ggsave and the people that created it; the export function on Rstudio works
# perfectly fine, ggsave is being a pain...

#ggsave("SAR-Flg22-1-2025.pdf", 
#       path = "", 
#       dpi = 300, 
#       width = 6, height = 8, units = "in")