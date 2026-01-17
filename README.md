#Update Week of Jan 11, 2026 -- Assignment 2

New scripts have been added:
     Carella-cleanup.R -------- Cleaning up the 2 Carella phloem proteomes and 
                                  saving as .rds files as per instructions
     hunting-for-anomalies.R -- Creating histograms using the mean abundances
                                  to look for anomalies, and identifying/analyz-
                                  ing any any extreme values in the data

<!-- BMB: what is RND? do you mean rds? -->

The .RND files are found in 'data/' as phloem_1_sum.rds, and phloem_2_sum.rds
The scripts are found in 'scripts'

I've used .gitignore to ignore a ~100 mb file I downloaded to potentially add 
GO features, or location features to the proteome(s) ('data/tair.gaf')

Some data has also been added including 'kvitko-proteome.csv' which the results
 from LC-MS/MS of the Arabidopsis apoplast following induction of a plant immune
 pathway called PAMP-triggered immunity. This is taken from Chen et al., 2025
 
Reference: Chen, H‑C., Newton, C. J., Diaz, G., Zheng, Y., Kong, F., Yao, Y.,
Yang, L., & Kvitko, B. H. (2025). Proteomic snapshot of pattern‑triggered 
immunity in the Arabidopsis leaf apoplast. The Plant Journal.

Goals/Objectives of this data analysis:
  The main goal is to determine similarity between proteins observed in the 
  apoplast of plants with an active immune response and compare this to 
  proteins that have been observed in the phloem of SAR-induced plants. The 
  premise of this idea is that proteinaceous components in the apoplast of the 
  infected leaf that may function in local immune signaling, may overlap with 
  proteinaceous components in the phloem that signal distant leaves to induce 
  resistance during SAR. I have found two apoplast phloem proteomes which I can 
  make comparisons to. 
   
  A separate goal I have is to potentially make reasonable/defensible calls on 
  proteins that are significantly abundant in the SAR proteomes but have been 
  disincluded based on the 2 peptide rule. I have briefly read on the use of the 
  2 peptide rule in proteomics data, and believe it may be possible to include
  proteins called based on single peptide matches, based on the quality of the 
  call, the size of the protein, and how consistently the peptide shows up across
  replicates. This is something I need to read more about, but may include some 
  interesting proteins that are lost based on filtering of proteins called with 
  only 1 unique peptide. 

## Week of Jan 4, 2026

New to this so not sure what the conventions are


#DIRECTORIES
data > Includes proteomics data from two independent experimental replicates, already normalized
       Data was collected by Dr. Phil Carella, formerly of the Cameron Lab (my lab)
       Published in Plant Physiology April 2016 (DOI: 10.1104/pp.16.00269)

scripts > Includes R scripts

Phloem_Proteome files:
   Accession -- Accession number for Arabidopsis protein
   Peptides used for quantitation -- Number of unique peptides used for quantitation
   Confidence -- Confidence score given by Progenesis QI software
   ANOVA (p) -- ANOVA test results between the groups (mock, vir, avr)
   Description -- I believe TAIR descriptions for the protein
   VLC.../DES... -- normalized abundance values for biological replicates 
         ends_with("MockPEX#") : protein abundance from plants that were infiltrated with mock inoculum (MgCl2)
         ends_with("VirPEX#") : protein abundance form plants that were infiltrated with Pst DC3000
         ends_with("AvrPEX#") : protein abundance from plants that were infiltrated with Pst expressing avrRpt2

#Context for the dataset:

The purpose of the proteomics was to identify proteins in the phloem exudates (i.e., phloem sap) of plants 
induced for a form of plant immunity  called Systemic Acquired Resistance (SAR). SAR is a immune response in 
plants where following infection of a leaf, long-distance signals are generated/activated to move systemically
via the phloem (vascular tissue) to distant leaves, where they signal the distant uninfected leaves to become
primed for defense, prior to interaction with any pathogens. This results in a more rapid and robust defense reponse
to subsequent pathogens. The phloem sap from mock-, vir-, or avr-induced plants was collected and subjected to liquid 
chromatography tandem mass-spectrometry to identify differentially abundant proteins in the phloem sap. 

#Biological Questions:

Initially I took this dataset as we had intended to perform additional LC MS/MS in our lab using different methods 
of SAR induction which did not rely on a live pathogen. We planned to induce SAR by expressing a protein that 
is recognized by plant immune receptors under and estrogen-inducible promoter, or infiltrating with a bacteria-derived 
peptide that is recognized by other plant immune receotors, to induce two distinct but overlapping immune pathways. We
then wanted to collect exudates from these plants for comparative proteomics. 

The question we wanted to ask was "Following SAR induction, what proteins enter the phloem in Arabidopsis"
More specifically we were asking "What proteins are more abundant in the phloem of plants induced for PAMP-triggered immunity" 
                  and..          "What proteins are more abundant in the phloem of plants induced for effector-triggered immunity"
Since both PAMP- and Effector triggered immunity are sufficient to induce SAR, and Pst avrRpt2 induces both immune pathways
simultaneously, we reasoned that common proteins found in the phloem proteomes of each induction method would comprise of 
important SAR long-distance signaling components. Moreover, since we are using a pathogen-free system, we could avoid 
any effects that live pathogens and their secreted effector proteins could have on the proteinaceous contents of the 
phloem during SAR signaling. 

Due to renovations and the building of the greenhouse, we could not perform the large-scale phloem sap collection
we intended since plants are incredibly sensitive, and stress can induce immune responses that would contaminate our 
results. Instead of abandoning the idea, I decided it might still be worth doing comparisons of our dataset to 
other publicly available SAR proteome datasets to see if any interesting proteins are revealed. I have found two
intercellular/apoplast SAR proteome datasets which I might like to compare our phloem proteome to. 

#Note: this might not be the dataset I continue with for the course, but it is one I could easily access and was familiar to me. 

Our lab also has an RNA-seq dataset which has not been well analyzed, and I may opt to use that to perform some analyses
once I have figured out a line of questioning pertinent to the data. I might want to compare the PTI dataset to published SAR
RNA seq data, RNA-seq data from plants induced in various ways. 
