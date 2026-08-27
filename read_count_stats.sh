#!/bin/bash

# Basecalling and Demultiplexing Script for Influenza DVG Drug Screen
# Script: demultiplexing.sh
# This file is shell script that downloads and processes raw sequencing reads
# Requires three arguments that should be passed from the command line:
#    run information, cross information, and experiment folders

#This script is part of the following manuscript:
#Influenza 
#Ilechukwu Agu  | Ivy José |  Samuel L. Díaz Muñoz

#Exporting Read Count Files
for file in *.mapped.sorted.BAM
do
	prefix=${file%.mapped.sorted.BAM}
	bam-readcount -b 20 -w 1 ${prefix}.mapped.sorted.BAM > ${prefix}.tab
done
