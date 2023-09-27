# Influenza DVG Drug Screen

This github repository includes code (and links to data) from the manuscript:  
Influenza A virus reassortment is strain dependent
Influenza DVG drug screen

If you are reading or using this, let us know how these data were useful for you. If you use these data and code, please cite the repository or the paper. Always open to collaborate! Please contact us!

Info below is templated, needs to be filled in 

### Quick Start
1. Make sure packages are installed (see #2 below) or use gbbseq-env.yml to set up Anaconda environment
2. git clone https://github.com/sociovirology/human_influenza_GbBSeq.git
3. chmod +x demultiplexing_GbBSeq6.sh
4. ./demultiplexing_GbBSeq6.sh runA "shared/cross_list_runA_pairwise.txt" pairwise_infections (pairwise)
5. ./demultiplexing_infection_conditions.sh (controls)
6. chmod +x amplicon_curation_strain_assignment.sh
6. ./amplicon_curation_strain_assignment.sh runA "shared/cross_list_runA_pairwise.txt" pairwise_infections
6. Rscript aiv_detection_environment_analysis.R (or load interactively in R)
7. Rscript aiv_detection_environment_analysis.R (or load interactively in R)

### CONTENTS
1. Project Description
2. Packages and software used to test code
3. Data
4. Code

### 1. Project Description
Influenza A viruses 

Abstract:
RNA viruses 

### 2. Packages and software used to test code
This code was tested using the following software packages:

* cutadapt (2.6)
* PEAR (0.9.6)
* usearch (8.1.1861)
* R (3.6.3 (2020-02-29) -- "Holding the Windsock") with packages:
    + dplyr, ggplot2, tidyr, reshape2, readr, gridExtra, ggthemes

Anaconda environment file is available in gbbseq-env.yml

### 3. Data
Data consists of sequencing output from the illumina MiSeq platform, sample information, reference database, exprimental coinfection titers, and sample barcodes

1) Sequencing file is available in the Sequence Read Archive (Accession SRX7014890)

2) Information on experimental coinfections is in data/cross_data_runA.csv 

3) Barcode information is in shared/barcodes3.fasta, shared/barcodes3rev.fasta shared/barcodes5.fasta, and shared/barcodes5rev.fasta

4) Database of reference sequences for the amplicons from influenza viruses used in the experimental coinfections is in shared/reference_database.fasta

5) Titers of experimental coinfection supernatants are in data/supernatant_titers.csv (Not needed to run code)

### 4. Code
Below are descriptions of the code files used to generate the tables, figures, and statistics in the paper.

1) demultiplexing_infection_conditions.sh: This file is shell script that downloads raw sequencing reads, demultiplexes each read; and generates a flat text summary file used in downstream analyses 

2) pairwise_infections.R: This file is an R script that analyzes reassortment in pairwise experimental coinfections among 5 human influenza A virus strains

3) control_infections.R: This file is an R script that analyzes reassortment in control experimental coinfections to test the methods and the GbBSeq appraoch