#Formal cleaning script; try to do a cleaner/better job at cleaning the data
# and think more about the data itself and how to treat it for use with downstream 
# processes. Try to write the code cleaner/more concise. 

#Libraries
library(tidyr)
library(stringr)
library(dplyr)

#Read in the first proteomics replicate
phloem_1_raw <- read.csv("data/phloem-proteome-1.csv", header = FALSE)

#Check the structure of the data.frame
head(phloem_1_raw)
str(phloem_1_raw)

#Things to address/clean:
#  [x]1. The headers are located on row 3, and there are header-headers(?) in rows 1-2
#  [x]2. The columns are all being read as chr despite being numerical
#  [x]3. There are rows with multiple Arabidopsis gene IDs (AT#G#####), need to 
#     address this somehow if I want to compare, or merge down the line. 
#  [x]4. There are calculated t-tests and ratios which I cannot trust; will need to
#     perform these in R and replace the rows
#  [x]5. The current headers are awful, regex will be annoying to figure out and they
#     currently contain machine-info, give them nicer descriptive names
#  [x]6. There shouldnt be any NA in any of the columns, but I should check to make sure
#      since it should be easy enough to do. Perhaps I can write a function to parse
#      each column of a data.frame for NA values, and report which columns contain 
#      NA values, maybe print the rows for quick reference.This would be nice to have
#      for cleaning other data.frames in the future. 


#___________Cleaning up the data: putting the data in frame_____________________

#Removing the first 2 rows by indexing the original dataframe
# Real dimensions of the data is 4:nrow(phloem_1_raw), 1:20

phloem_1 <- phloem_1_raw[4:nrow(phloem_1_raw), 1:20]

#Rename the headers to the contents of row 3 and rename the headers to something
# that is easier to read and type. Easy solution is to just create a character 
# vector with all the names in order, but I want to try and remove the "VLC..." 
# and keep the ends (i.e., MockPEX#) instead. (the R way)


col_names <- phloem_1_raw[3, 1:20] %>% 
  str_replace(pattern = ".*-", replacement = "")

#Update the column names to these...

colnames(phloem_1) <- col_names

#Check my work so far
colnames(phloem_1)
str(phloem_1)
head(phloem_1)

#Looking much better; need to fix the column types now; set all the appropriate 
# columns to be numeric rather than chr. Columns to be changed: Peptides used for Quant., 
# Confidence Score,  Anova (p), Mock/Vir/AvrPEX#, ratio/ttest.

#This is the third iteration of this process I have tried; much cleaner/compact code

phloem_1 <- phloem_1 %>% 
  mutate(across(c(`Peptides used for quantitation`, 
                  `Anova (p)`, 
                  `Confidence score`, 
                  contains("PEX")), 
                ~as.numeric(.x)))

str(phloem_1)

#Check if any of my columns contain any NA values... 

#Try to create a function to quickly check each column for NA values. Try to 
# create a discrete function to do this

AnyNA_Columns <- function(data) {
  for (i in names(data)) {         # For i in names(data) (iterate through columns)
    if (any(is.na(data[i]))) {     # If there are any NA values in data[i], 
      message(i, " HAS NA VALUES!??!?!??!!!") # print to console "column name" "has NA values"
    } else {                       # else 
      message(i, " CLEAR")         # print "column name" "clear"
    }
  }
}

AnyNA_Columns(phloem_1)

#Create new columns to replace old calculations and summarize the important data

phloem_1_sum <- phloem_1 %>% 
  rowwise %>% 
  mutate(mean_mock  = mean(c_across(starts_with("MockPEX")), na.rm = TRUE), 
         mean_vir = mean(c_across(starts_with("VirPEX")), na.rm = TRUE), 
         mean_avr = mean(c_across(starts_with("AvrPEX")), na.rm = TRUE), 
         FC_vir = mean_vir/mean_mock, 
         FC_avr = mean_avr/mean_mock, 
         ttest_mock_vir = t.test(c_across(starts_with("MockPEX")), 
                                c_across(starts_with("VirPEX")), 
                                alternative = "less")$p.value, 
         ttest_mock_avr = t.test(c_across(starts_with("MockPEX")), 
                                c_across(starts_with("AvrPEX")), 
                                alternative = "less")$p.value,) %>% 
  select(Accession, `Peptides used for quantitation`, Description, 
         mean_mock, mean_vir, mean_avr, 
         FC_vir, FC_avr, ttest_mock_vir, ttest_mock_avr) %>% 
  ungroup()


str(phloem_1_sum)

#Dealing with Multiple gene IDs and the description column. 

#Data is now cleaned except for the fact there are some rows with multiple 
# accesssion identifiers. I'm not really sure the most appropriate way to 
# address this. Biologically they are essentially the same protein or at least 
# equally-likely proteins identified based on the peptides. I believe this tends
# to occur when there are highly similar proteins like homologs/paralogs arisen 
# from duplications. This appears to be quite standard practice. 

#Cleaning up the description column:
#  Description includes a bunch of variables which could probably be more easily 
#  Handeled by seperating them into seperate columns. 

phloem_1_sum[["Description"]][1:5]

#Contains the gene symbol, the name, the position
#Also appears to be some GO stuff; I think if I want to include any of this kind of 
# data in the data.frame, I'll just download the descriptions and gene symbols  
# then create a new column myself with this information. 

#Found two datasets from The Arabidopsis Information Resource (TAIR) 
#  One contains functional and curated descriptions for each protein
#  One contains gene symbols and also interpro IDs which might be useful
#  One of them is >50 mb, the other is only ~20 mb

tair_functional_desc <- read.table("data/TAIR10_functional_descriptions.txt", 
                                   sep = "\t",
                                   quote = "",
                                   fill = TRUE,
                                   header = TRUE)

#Two things I want to do with this dataframe before merging;
#  1. Remove the Type column 
#  2. Extract the first element of the string from the computational description 
#     which has the gene name and symbol

tair_functional_desc <- tair_functional_desc %>% 
  mutate(gene_name = str_extract(tair_functional_desc$Computational_description, "^[^;]*")) %>% 
  select(!c(Type, Computational_description))

#This is also where having the multiple accessions is going to bite me, try a quick
# and dirty solution by just cloning the rows with each accession ID, they will all
# show up together in any analysis I do so will be hard to assume they are seperate things...
# to be safe, I think I will add a new row called proteome_1_ID so that when I split them, 
# I can easily group them and treat them as a single thing by using the proteome_1_ID to index
# them. 

phloem_1_sum <- phloem_1_sum %>% 
  mutate(proteome_1_ID = 1:nrow(phloem_1_sum)) %>% 
  separate_rows(Accession, sep = ";")

#Merging dataframes by Accession

phloem_1_sum_desc <- left_join(x = phloem_1_sum,
                               y = tair_functional_desc,
                               by = c("Accession" = "Model_name")) %>% 
  select(!Description)

saveRDS(phloem_1_sum_desc, "data/phloem_1_sum.rds")

#Happy with the state of the data.frame, so will save the dataframe as a .csv
# so when I do comparisons I can begin from a clean data.frame in a fresh script

#write.csv(phloem_1_sum_desc, "data/phloem_1_sum.csv", col.names = TRUE, row.names = FALSE)


###########################Phloem_Proteome_2#####################################

#Cleanup environment before starting:

rm(phloem_1, phloem_1_a, phloem_1_c, phloem_1_raw, phloem_1_sum, phloem_1_sum_desc)

#Read in the second proteomics replicate... 
phloem_2_raw <- read.csv("data/phloem-proteome-2.csv", header = FALSE)


#Check structure and shape of the data.frame
str(phloem_2_raw)
head(phloem_2_raw)

#Headers are on row 3, first row of data is row 4

phloem_2 <- phloem_2_raw[4:nrow(phloem_2_raw), 1:ncol(phloem_2_raw)]

#Take row 3 as headers, and adjust header names for the abundance values

col_names <- phloem_2_raw[3,] %>% 
  str_replace("^[^-]*-", "")

#Give columns headers using col_names

colnames(phloem_2) <- col_names

#Fix column types:

phloem_2 <- phloem_2 %>% 
  mutate(across(c(`Unique peptides`, 
                  `Confidence score`, 
                  contains("PEX")), 
                ~as.numeric(.x)))

#Apparently NAs introduced by coercion

AnyNA_Columns(phloem_2)

#Those columns are going to be removed so it's alright to leave alone


#Perform calculations for means, ttests, summarize dataframe

phloem_2_sum <- phloem_2 %>% 
  rowwise %>% 
  mutate(mean_mock  = mean(c_across(starts_with("Mock-PEX")), na.rm = TRUE), 
         mean_vir = mean(c_across(starts_with("Vir-PEX")), na.rm = TRUE), 
         mean_avr = mean(c_across(starts_with("Avr-PEX")), na.rm = TRUE), 
         FC_vir = mean_vir/mean_mock, 
         FC_avr = mean_avr/mean_mock, 
         ttest_mock_vir = t.test(c_across(starts_with("Mock-PEX")), 
                                 c_across(starts_with("Vir-PEX")), 
                                 alternative = "less")$p.value, 
         ttest_mock_avr = t.test(c_across(starts_with("Mock-PEX")), 
                                 c_across(starts_with("Avr-PEX")), 
                                 alternative = "less")$p.value,) %>% 
  select(Accession, `Unique peptides`, Description, 
         mean_mock, mean_vir, mean_avr, 
         FC_vir, FC_avr, ttest_mock_vir, ttest_mock_avr) %>% 
  ungroup()

#Create index for proteome 2 and seperate rows with multiple gene_ids
phloem_2_sum <- phloem_2_sum %>% 
  mutate(proteome_2_ID = 1:nrow(phloem_2_sum)) %>% 
  separate_rows(Accession, sep = ";")

#Merging with tair_functional descriptions and removing the description column
phloem_2_sum_desc <- left_join(x = phloem_2_sum,
                               y = tair_functional_desc,
                               by = c("Accession" = "Model_name")) %>% 
  select(!Description)


#saveRDS(phloem_2_sum_desc, "data/phloem_2_sum.rds")

#Data.frame looks good; save cleaned dataframe
#write.csv(phloem_2_sum_desc, "data/phloem-2-sum.csv", col.names = TRUE)
