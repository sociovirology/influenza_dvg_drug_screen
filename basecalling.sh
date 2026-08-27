#!/bin/bash

# Basecalling Script for Influenza DelVG Drug Screen
# Script: basecalling.sh
# This file is shell script that takes the ONT (Oxford Nanopore Technologies) sequencing run folder output and basecalls FAST5 files to generate FASTQ files that downloads and processes raw sequencing reads
 #In practive this file is not run by users as the FAST5 files are too large to be provided online.
 #This file is provided here to provide insight into how the sequence data was generated (pre-FASTQ)

# Requires entire ONT output directory.

#This script is part of the following manuscript:
#Influenza A deletion-containing viral genome production is altered by metabolites, metabolic signaling molecules, and cyanobacterial extracts 
#Ilechukwu Agu | Suzanna Y. Gomez | Ivy R. José | Samuel L. Díaz-Muñoz


#Basecall command
~/ont-guppy/bin/guppy_basecaller --input_path '/media/user/Lab Portab/dvg_dug_screen_ile_v2/no_sample/20230921_1151_MN23913_FAW91824_5830ce8e/' --save_path /mnt/data0/sam/dvg_drug_screen_ile/fastq_hac/ -c dna_r10.4.1_e8.2_400bps_hac.cfg -x cuda:all:100% --chunks_per_runner 1024 --recursive --calib_detect --do_read_splitting

#ONT Guppy basecalling software version 6.5.7+ca6d6af, minimap2 version 2.24-r1122
#config file:        /home/user/ont-guppy/data/dna_r10.4.1_e8.2_400bps_hac.cfg
#model file:         /home/user/ont-guppy/data/template_r10.4.1_e8.2_400bps_hac.jsn
#input path:         /media/user/Lab Portab/dvg_dug_screen_ile_v2/no_sample/20230921_1151_MN23913_FAW91824_5830ce8e/
#save path:          /mnt/data0/sam/dvg_drug_screen_ile/fastq_hac/
#chunk size:         2000
#chunks per runner:  1024
#minimum qscore:     9
#records per file:   4000
#num basecallers:    4
#gpu device:         cuda:all:100%
#kernel path:        
#runners per device: 4

#Use of this software is permitted solely under the terms of the end user license agreement (EULA).
#By running, copying or accessing this software, you are demonstrating your acceptance of the EULA.
#The EULA may be found in /home/user/ont-guppy/bin
#loading new index: /home/user/ont-guppy/data/lambda_3.6kb.fasta
#Found 21294 input read files to process.
#Init time: 24981 ms
#
#0%   10   20   30   40   50   60   70   80   90   100%
#|----|----|----|----|----|----|----|----|----|----|
#***************************************************
#Caller time: 24913950 ms, Samples called: 263620048773, samples/s: 1.05812e+07
#Finishing up any open output files.
#Basecalling completed successfully.
```
#Results
#FASTQ HAC Fail
#26485808 total lines
#*6,621,452 reads*
#FASTQ HAC Pass
#62210552 total lines
#*15,552,638 reads*
#*Pass read percentage 70.14% out of 22,174,090 total reads *

#Let us know it's done
echo "Finished!"