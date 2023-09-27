#!/bin/bash

# Basecalling and Demultiplexing Script for Influenza DVG Drug Screen
# Script: demultiplexing.sh
# This file is shell script that downloads and processes raw sequencing reads
# Requires three arguments that should be passed from the command line:
#    run information, cross information, and experiment folders

#This script is part of the following manuscript:
#Influenza 
#Ilechukwu Agu  | Ivy José |  Samuel L. Díaz Muñoz

#Verify we have 3 arguments

if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments: run, cross list, and experiment"
    echo "run has format 'runA'"
    echo "File name for cross list. Cross list file is just a plain text file one line per cross: cross1 [break] cross9 etc. "
    echo "Experiment can be any folder name you desire, e.g. 'infection_conditions'"
    exit
fi

#Checking command line parameters:
echo "Checking command line arguments: run, cross list, and experiment"
echo $1
echo $2
echo $3

#Download FASTQ sequencing files from SRA via EBI (Illumnia PE read files)
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/043/SRR21880043/SRR21880043_1.fastq.gz -O data/runA_R1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/043/SRR21880043/SRR21880043_2.fastq.gz -O data/runA_R2.fastq.gz

#Unzip files
#gunzip data/runA_*.fastq.gz

#First let's basecall
~/ont-guppy/bin/guppy_basecaller --input_path /mnt/data0/sam/dvg_drug_screen_ile/ --save_path /mnt/data0/sam/dvg_drug_screen_ile/quick_check/fastq_hac --config dna_r10.4.1_e8.2_400bps_hac -x cuda:all:100% --chunks_per_runner 1024 --recursive --calib_detect --do_read_splitting

#Then let's demultiplex

#Let us know it's done
echo "Finished!"