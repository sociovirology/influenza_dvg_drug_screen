#!/bin/bash

# Demultiplexing Script for Influenza DelVG Drug Screen
# Script: demultiplexing.sh
# This file is shell script that takes FASTQ's basecalled by Guppy (see basecalling.sh )  ONT (Oxford Nanopore Technologies) sequencing run folder output and basecalls FAST5 files to generate FASTQ files that downloads and processes raw sequencing reads
 #In practive this file is not run by users as the demultiplexed FASTQ's are downloadable from ENA/SRA (SEE FILENAME)
 #This file is provided here to provide insight into how the demultiplexed FASTQ's were generated

# Requires entire ONT output directory and basecalling.sh to be run first.

#This script is part of the following manuscript:
#Influenza A deletion-containing viral genome production is altered by metabolites, metabolic signaling molecules, and cyanobacterial extracts 
#Ilechukwu Agu | Suzanna Y. Gomez | Ivy R. José | Samuel L. Díaz-Muñoz


#Demultiplexing commands 

#Require both barcodes minimum scores F: 60, R: 40 

~/ont-guppy/bin/guppy_barcoder --require_barcodes_both_ends --input_path /mnt/data0/sam/dvg_drug_screen_ile/fastq_hac/pass --save_path /mnt/data0/sam/dvg_drug_screen_ile/fastq_hac/pass/double_barcodes_60_40/ --detect_mid_strand_barcodes --barcode_kits EXP-PBC096 --min_score_barcode_rear 40 --min_score_barcode_front 60 -x cuda:all:100% --worker_threads 32
#ONT Guppy barcoding software version 6.5.7+ca6d6af
#input path:         /mnt/data0/sam/dvg_drug_screen_ile/fastq_hac/pass
#save path:          /mnt/data0/sam/dvg_drug_screen_ile/fastq_hac/pass/double_barcodes_60_40/
#min. score front:   60
#min. score rear:    40
#gpu device:         cuda:all:100%

#Found 5324 input files.

#0%   10   20   30   40   50   60   70   80   90   100%
#|----|----|----|----|----|----|----|----|----|----|
#***************************************************
#Done in 18009466 ms.

#Let us know it's done
echo "Finished!"