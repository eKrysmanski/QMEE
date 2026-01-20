#Plots to check data for anomolies 

#Libraries
library(ggplot2)
library(patchwork)

#reading in the data
phloem_1 <- readRDS("data/phloem_1_sum.rds")
phloem_2 <- readRDS("data/phloem_2_sum.rds")


#Quick question for Dr. Dushoff/Bolker: I tend to make a lot of objects in my
# scripts, but it can obviously fill up the environment with things that may be 
# temporary or just to check my work. Is this just a normal part of developing a 
# script which is later removed, or is there a better/cleaner way to write.

## BMB: making lots of temporary objects is absolutely fine.
## When you chunk your workflow into separate pieces and restart R/
## run each piece as a separate R processes (and set RStudio up to
## *not* automatically save/restore the global environment (workspace),
## all the temporary objects will be transient ...

#Let's check log-distribution of mean protein abundances; should 
# be log-normal, and hopefully reveal any strange values

scatter_p1_mock <- ggplot(data = phloem_1, aes(x = mean_mock + 1)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,60))

scatter_p1_vir <- ggplot(data = phloem_1, aes(x = mean_vir + 1)) + 
  geom_histogram(bins = 50, fill = "green") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,60))

scatter_p1_avr <- ggplot(data = phloem_1, aes(x = mean_avr + 1)) +
  geom_histogram(bins = 50, fill = "orange") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,60))

#patchwork lets me facet; show three at the same time side-by-side
(scatter_p1_mock | scatter_p1_vir | scatter_p1_avr)

## BMB: you can generate variants of plots like this, without
##  repeating as much code (although this doesn't override the colour)
scatter_p1_vir <- scatter_p1_mock + aes(x = mean_vir + 1)
scatter_p1_avr <- scatter_p1_mock + aes(x = mean_avr + 1)

## BMB: however, the idiomatic way to do this is to pivot to longer:
phloem_1_long <- phloem_1 |> pivot_longer(cols = starts_with("mean"))

gg_facet <- ggplot(data = phloem_1_long, aes(x = value+1, fill = name)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  facet_wrap(~ name) +
  scale_fill_manual(values = c("steelblue", "green", "orange"))

print(gg_facet)

## BMB: note, fill colour is pretty but redundant ...

#There are infinite values introduced by log10; this can occur if there are 
# values that are 0, or values that are infinity. If there were something
# strange like a negative value then there would be an error. Fix by adding a 
# "pseudocount" of + 1 to all values.

## BMB: this is potentially a much deeper question - probably OK for
## exploration, and probably OK if the values are integers

#This all looks fine; proteomics data is supposed to look log-normal...

## BMB: DRY. You can "add" a new data set to an existing ggplot ...
phloem_2_long <- phloem_2 |> pivot_longer(cols = starts_with("mean"))
gg_facet + phloem_2_long

#Check out proteome-2
scatter_p2_mock <- ggplot(data = phloem_2, aes(x = mean_mock + 1)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,100))

scatter_p2_vir <- ggplot(data = phloem_2, aes(x = mean_vir + 1)) + 
  geom_histogram(bins = 50, fill = "green") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,100))

scatter_p2_avr <- ggplot(data = phloem_2, aes(x = mean_avr + 1)) +
  geom_histogram(bins = 50, fill = "orange") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,100))

(scatter_p2_mock | scatter_p2_vir | scatter_p2_avr)

#This also appears reasonable


############################################################################
#######         Checking very large and very small values       ############

#This also looks reasonable

#There appears to be a couple of very small and very large values for each
# let's extract those and see why they are so extreme

## BMB: again, easier to do with long-format data
phloem_1_low <- phloem_1 %>% 
  filter(mean_mock < 100 | mean_vir < 100 | mean_avr < 100)

phloem_1_low

#Nothing strange here

phloem_1_high <- phloem_1 %>% 
  filter(mean_mock > 5E7 | mean_vir > 1E7 | mean_avr > 5E7)

phloem_1_high

#Nothing strange here, these are high across all treatments and not significant

## DRY: write a checking function?

#Do the same for phloem_2
phloem_2_low <- phloem_2 %>% 
  filter(mean_mock < 100 | mean_vir < 100 | mean_avr < 100)

phloem_2_low

#A couple of the FC are inf, since mock shows "0" abundance, but there are 
# proteins in the vir and/or avr. 

phloem_2_high <- phloem_2 %>% 
  filter(mean_mock > 5E7 | mean_vir > 1E7 | mean_avr > 5E7)

phloem_2_high

#Nothing strange here

######          Finding any infinite values in the data.frames    ###############

## BMB: again, easier in long format
#For phloem_1...
p1_inf <- phloem_1 %>%
  filter(if_any(where(is.numeric), is.infinite))

p1_inf

#For phloem_2...
p2_inf <- phloem_2 %>% 
  filter(if_any(where(is.numeric), is.infinite))

p2_inf

#Quick check, do they show up in both proteomes?

#Determine the IDs
inf_proteins <- as.vector(p2_inf$Accession)

#Check for these in phloem_1

phloem_1_inf <- phloem_1 %>% 
  select(Accession == any_of(inf_proteins))

phloem_1_inf <- phloem_1 %>%
  filter(Accession %in% inf_proteins)

                 
#####Notes:                                        
# So these infinite values show up only in phloem_2...
# How do I deal with these, and do I need to deal with these:
#     On one hand, it is standard to add a pseudocount (+1) and log2() transform
#     On the other hand, biologically this shows that these proteins are present in 
#     the vir/avr proteomes, but not the mock and are therefore potentially quite significant

# I think I will leave it for now, but it might bite me in the butt later on

# mark: 2.1
