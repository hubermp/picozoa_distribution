_______________________________________________________________

### Global Distribution, Diversity, and Ecological Niche of Picozoa, a Widespread and Enigmatic Marine Protist Lineage

This repository contains the data files used to explore the distribution and diversity of Picozoa in the EukBank database. The data is organized into several files, each serving a specific purpose in the analysis of the global distribution, diversity, and ecological niche of Picozoa.

#### Main R Script
Picozoa_distribution_and_diversity.Rmd

_______________________________________________________________

#### Data Description
##### Data Used to Study the Global Distribution and Diversity of Picozoa

1. EukRibo_taxonomic_EukBank_db.csv  
This table contains the relative contribution of high-rank protistan groups to the total number of reads in the EukBank dataset, computed using the EukRibo reference database.
2. metadata_associated_TableS1.csv  
This table contains the contextual data of the samples from the EukBank dataset analyzed in this work. Please refer to Supplementary Table S1 in the manuscript.
3. abundance_table_marine_samples-rarefy.csv  
This table contains the number of reads (post-rarefaction) for each pOTU in EukBank pico-size marine samples, excluding time series and replicate data.

##### Data Used to Study the Phylogenetic Community Structure and Phylogenetic Niche Conservatism of Picozoa

4. phylo_abundance_table_marine_samples-rarefy.csv  
This table contains the number of reads (post-rarefaction) for each pOTU in EukBank pico-size marine samples, excluding time series and replicate data, along with reference sequence IDs used to construct the phylogenetic tree, with an abundance of 0 used for MNTD analysis.
5. environmetal_variables_WOA.csv  
This table contains the environmental variables obtained from WOA, used for determining the pOTU environmental optima.
6. RAxML_bipartitions.Picozoan.tree
This file contains the Picozoa 18S rDNA maximum likelihood tree. The tree was constructed using reference sequences and pOTUs, considering 1 000 replicates. 
7. PhyloDist_NichOverlap.csv
This table contains the phylogenetic distance and niche overlap (Schoener’s D metric) values between each abundant pOTU.
