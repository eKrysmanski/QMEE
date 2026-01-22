#### Barplot for SAR Assays -- Full Labs + Stats
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(showtext)

#Read in data
data <- read.csv("data/SAR-assay-data/SAR-flg22-1-2025-ANOVA.csv")

################################################################################
#####################      Setting up factors       ############################

#Specify columns as factors with levels as necessary
#Note that this experiment only had Col-0 plants, but others may have mutants 
# which is why it is being included here. 
data$treatment <- factor(data$treatment, levels = c("MgCl2", 
                                                    "Pst AvrRpt2",
                                                     "H2O",
                                                    "Flg22"))

data$induction <- factor(data$induction, levels = c("Mock", "SAR"))

data$athal_genotype <- factor(data$athal_genotype, levels = c("Col-0"))

data$group <- factor(data$group, levels = c("1", "2"))

################################################################################
##################### Misc. bits and pieces for plot ###########################

#Load in google font--------------------------------------------------------
font_add_google(name = "Pt Serif", family = "pt_serif")  #Adding google font
showtext_auto()

#Create a summarized data.frame because I find it easier for adding jitterpoints
# when I can use the summarized data for the major geoms

data_sum <- data %>% 
  group_by(athal_genotype, treatment, group, induction) %>% 
  summarize(mean_pop = mean(population), 
            sd_pop = sd(population)) %>% 
  ungroup()

################################################################################
###  Inputs: Can be adjusted to quickly modify major components of the plot  ###
################################################################################

#y-axis inputs--------------------------------
y_lims <- c(1e5, 5e7)
y_breaks <- 10^(5:7)

#Labels----------------------------------------
bar_labs <- c(
  expression(MgCl[2]),                  #Note: to stack use "expression(atop("text", MgCl[2]))
  expression(italic("Pst AvrRpt2")),
  expression(H[2]*O),
  expression(Flg22)
)

#group labels-----------------------------------
group_labs <- c("Col-0", "Col-0",  #This is a bit broken, for some reason needs a vector length 4 
                "Col-0", "Col-0")  #Works if you just specify each label identically twice since they stack

#Stats------------------------------------------
stats_labs <-c("a", "c", "a", "b") #One-way ANOVA + Tukey's HSD

################################################################################
############################ Plot Script #######################################

plot <- ggplot(data_sum, aes(x = group, y = mean_pop, fill = induction)) +
  geom_col(                                                    #Columns
    position = position_dodge(0.8),                              #Dodge
    width = 0.7,                                                 #Essentially bar size/gap within groups
    colour = "black",                                            #Border colour
    linewidth = .9                                               #Border width
    ) +                                            
  geom_errorbar(                                               #Adding standard deviation bars
    position = position_dodge(0.8),                              #Need to dodge
    aes(x = group,                                               #Aesthetic mappings
        ymax = mean_pop + sd_pop,                                #Positive whiskers
        ymin = mean_pop - sd_pop),                               #Negative whiskers
    width = 0.2                                                  #width of the hats
    ) + 
  geom_jitter(                                                 #Adding jitterpoints
    data = data,                                                 #using raw dataframe
    position = position_jitterdodge(dodge.width = 0.8),          #special position_dodge for jitterpoints
    aes(x = group, y = population, fill = induction),            #Same aesthetics as rest
    alpha = 0.75                                                 #Transparency in case they overlap
    ) +
  scale_y_log10(                                               #Adding a logarithmic scale
    limits = y_lims,                                             #Specify axis limits (set above)
    breaks = y_breaks,                                           #Specify label breaks (set above)
    labels = trans_format("log10", math_format(10^.x)),          #Adds superscript for labels
    oob = scales::oob_squish,                                    #Prevents out-of-bounds data from being removed
    expand =  expansion(0)                                       #Prevents cushion between bottom of bars and axis
    ) +
  geom_text(             #BAR LABELS                           #Adding bar-lables
    aes(x = group,                                               #map to groups
        y = 1,                                                   #Y-coordinate for label (Set to 1 so I can vjust below axis)
        label = bar_labs),                                       #The labels (set above)
    family = "pt_serif",                                         #Font
    size = 6,                                                    #Font size 
    position = position_dodge(0.8),                              #Need to dodge
    vjust = 1.5                                                  #vertical adjust location of text
  ) +
  geom_text(             #GROUP LABELS                         #Group labels
    aes(x = group,                                                 #Same as BAR LABELS
        y = 1,                                                     #No position_dodge()
        label = group_labs),
    family = "pt_serif", 
    size = 6, 
    vjust = 4.5                                                    #Adjust even lower
  ) +
  geom_text(            #Label Bar                               #I acknowledge this is proabably very stupid, and
    aes(x = group,                                                  #there is probably a better way to get  line segments...
        y = 1,                                                      #but it works and is very easy to adjust as needed.  
        label = "_______________",                                  #I have tried geom_hline(), annotate(), geom_segment()
        fontface = 2),
    family = "pt_serif", 
    size = 7, 
    vjust = 1.95
  ) + 
  geom_text(            #STAT LABELS                  #Adding in significance letters from Anova + Tukey's HSD    
    aes(x = group,                                    
        y = mean_pop + sd_pop * 1.5,                    #Set y-coordinate at a fixed distance above each positive whisker
        label = stats_labs),                            #Labels are specified above
    family = "pt_serif", 
    size = 7, 
    position = position_dodge(0.8),
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
  scale_fill_manual(                                           #Create discrete scale (i.e. bar colours)
    values = c("Mock" = "#35B749","SAR" = "lightgoldenrodyellow")          #Think I like these colours;
    ) +                  #349433          #FFC107
  labs(
    y = "Bacterial density in distant leaves (cfu/ld)"
  ) +
  theme_minimal(base_size = 14) +                  #Basic theme default
  theme(                                           #Adjust the theme however you like
    text = element_text(family = "pt_serif"),
    panel.grid = element_blank(),                         #Remove grid lines
    axis.line = element_line(color = "black",             #Setting solid line for axes
                             linewidth = 0.9),              #Setting line width for axes
    axis.line.x = element_line(color = "black", 
                               linewidth = 0.9),
    axis.text.y = element_text(vjust = 0.2,               #Adjusting position and colour of y-labs
                               hjust = 2,                   #Horizontal adjustments
                               color = "black",             #Font colour
                               size = 17),                  #Font size
    axis.title.y = element_text(margin = margin(r = 15),  #Increase right margin
                                size = 20),                 #Font size
    axis.title.x = element_blank(),                       #Remove the x-axis title,
    axis.text.x = element_blank(),                        #Remove x axis text
    plot.margin = margin(t = 20, r = 30, b = 50, l = 40), #Set plot margins
    legend.position = "none"                              #Remove legend
  )

plot

#Saving Plot

ggsave("SAR-Flg22-1-2025.pdf", 
       path = "C:/Users/eKrys/Desktop/Cameron Lab/Cameron Lab - June 2025/2025 Experiments/SAR-SA-PEX-2026/SAR-SA-PEX-1-2026/", 
       dpi = 300, 
       width = 6, height = 8, units = "in")

################################################################################
############################  Playing with colours #############################

#Note: make sure whatever is used is colourblind accessible 

#Colour combos:
"#00A08A"  "#F2AD00"
"#00A08A"  "#F2AD00"
"#00A08A"  "#F2AD00"

"springgreen2"  "yellow2"
"springgreen2"  "yellow2"
"springgreen2"  "yellow2"

"#43ff43"  "#ff43a1"
"#43ff43"  "#ff43a1"
"#43ff43"  "#ff43a1"

"#66ffcc"  "#ffcc66"
"#66ffcc"  "#ffcc66"
"#66ffcc"  "#ffcc66"

"#3eff03"  "#ff033e"  "#03ffc4"
"#3eff03"  "#ff033e"  "#03ffc4"
"#3eff03"  "#ff033e"  "#03ffc4"

"#008000" "#800000"
"#008000" "#800000"
"#008000" "#800000"

"#12ffa8" "#ffa812" "#ffe312"
"#12ffa8" "#ffa812" "#ffe312"
"#12ffa8" "#ffa812" "#ffe312"

"#b2ec5d" "#eeee90"
"#b2ec5d" "#eeee90"
"#b2ec5d" "#eeee90"

"#78C841" "#FF9B2F"
"#78C841" "#FF9B2F"
"#78C841" "#FF9B2F"

"#349433" "#FFC107"
"#349433" "#FFC107"
"#349433" "#FFC107"

"#35B779" "#FDE725" "#31688E"
"#35B779" "#FDE725" "#31688E"
"#35B779" "#FDE725" "#31688E"

#From RBG Chrysanthemum 
"#627338" "#fe9a0f"
"#627338" "#fe9a0f"
"#627338" "#fe9a0f"

#From RBG Rose
"#78C841"   "#fea55f"
"#78C841"   "#fea55f"
"#78C841"   "#fea55f"

#From chlorotic arabidopsis 
"#acd800"  "#feef00"
"#acd800"  "#feef00"
"#acd800"  "#feef00"

#Deciding on Colours
Greens: 
  "green", "darkgreen", "forestgreen", "lightgreen", 
"mediumseagreen", "springgreen3", "yellowgreen", "darkolivegreen4"

Yellows: 
  "yellow", "lightyellow", "yellow1", "yellow2", "yellow3", 
"yellow4", "greenyellow", "lightgoldenrodyellow"

show_col("greens")
show_col("oranges")

#wesanderson palettes
library(wesanderson)

names(wes_palettes)
show_col(wes_palette("Darjeeling1"))


