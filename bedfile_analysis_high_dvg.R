#Recombination Plots and junctions
#CA09 samples only
#
#
#

#Load libraries
library(circlize)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(GGally)


#Now recombination plot start from scratch to include only select barcodes
#Get the list of bedfiles
file_names <- list.files("/Users/mixtup/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/output/adjusted_bedfiles_no_readnames_may2024", "*_Virus_Recombination_Results.bed", full.names = T)
file_names <- file_names[-7]

recombinations_bedfile <- NULL
for (file in file_names) {
  print(file)
  df <- read.table(file, skip=1, sep="\t", stringsAsFactors=FALSE, quote="")
  df <- cbind(df, rep(file, nrow(df)))
  recombinations_bedfile <- rbind(recombinations_bedfile, df)
}

#Name columns
colnames(recombinations_bedfile) <- c("reference", "start", "stop", "event_type", "count", "sense", "5_reads", "3_reads", "5_ref_sequence", "3_ref_sequence", "sample")

#Clean up the filename, remove path
recombinations_bedfile$sample <- gsub("/Users/mixtup/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/output/adjusted_bedfiles_no_readnames_may2024/","", recombinations_bedfile$sample)
recombinations_bedfile$sample <- gsub("_1-25_adjusted_Virus_Recombination_Results.bed","", recombinations_bedfile$sample)

#Note this doesn't include isolates, so need to update code to compare across isolates
recombinations <- recombinations_bedfile %>% group_by(reference, start, stop, event_type) %>% 
  mutate(number_of_samples=n()) %>% 
  mutate(count=sum(count %>% as.numeric())) %>% 
  ungroup()

#Also include a size column. #Calculating size of deletions. Taking absolute value of start-stop to account for negative strand events
recombinations <- mutate(recombinations, size = abs(stop-start))

#How many non-deletions
nrow(recombinations) - nrow(subset(recombinations, event_type == "Deletion"))
#[1] 110
#Keep only deletions
recombinations <- subset(recombinations, event_type == "Deletion")

#Plot
ggplot(recombinations, aes(x = stop, y = start, color = number_of_samples)) +
  geom_point(aes(size = count),alpha = 0.6) +
  scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference)


#What I really want is to see the difference in DVGs between samples / treatments

sample_data <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/data/sample_data.csv")
demultiplexing_trimming_stats <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/output/demultiplexing_trimming_stats.csv")

#Adjust columns to match for join
colnames(sample_data) <-  c("sample","sample_name", "well_position", "control", "group", "strain_name", "strain", "treatment", "bioreplicate", "cdna_concentration")

#Join sample
recombinations <- left_join(recombinations, sample_data, by = c("sample"))

#Junctions per sample
ggplot(recombinations, aes(x = treatment, color = sample_name)) +
  geom_bar() +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference)

#Junctions per sample
ggplot(recombinations, aes(x = reference, color = sample_name)) +
  geom_bar() +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ treatment)

#Junctions per sample with only >1 RSC
ggplot(subset(recombinations, count > 1), aes(x = treatment, color = sample_name)) +
  geom_bar() +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference)


#Junctions per sample, only cannonical junctions
ggplot(subset(recombinations, size > 190), aes(x = treatment, color = sample_name)) +
  geom_bar() +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference)
#Alright *there* we are. Huge difference in cannonical junctions with PA *in the polymerase segments*

#Does that change dramatically if I adjust RSC?
ggplot(subset(recombinations, size > 190 & count > 1), aes(x = treatment, color = sample_name)) +
  geom_bar() +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference)
#Not really! Looks quite good in relative terms. 


#Let's do a map of deletions like in Alnaji et al. Parallel coordinates plot
ggplot(recombinations, aes(x = treatment, color = sample_name)) +
  geom_bar() +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_grid(reference ~ treatment)

ggplot(subset(recombinations, reference == "PB2|Segment:1"), aes(x = stop, y = start, color = number_of_samples)) +
  geom_point(aes(size = count),alpha = 0.6) + 
  scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_grid(reference ~ treatment)

ggplot(recombinations, aes(x = stop, y = start, color = number_of_samples)) +
  geom_point(aes(size = count),alpha = 0.6) + 
  scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_grid(reference ~ treatment)

#By sample.
ggparcoord(subset(recombinations, reference == "PB2|Segment:1"), columns = 2:3, scale = "globalminmax") + facet_wrap(~ treatment)
ggparcoord(subset(recombinations, reference == "PB1|Segment:2"), columns = 2:3, scale = "globalminmax") + facet_wrap(~ treatment)

#Isolate to one segment and one case
ggparcoord(subset(recombinations, size > 190 & reference == "PB2|Segment:1" & sample == "barcode56"), columns = 2:3, scale = "globalminmax")
ggparcoord(subset(recombinations, size > 190 & reference == "PB1|Segment:2" & sample == "barcode56"), columns = 2:3, scale = "globalminmax")

#Issue here is that I need to flip the start and stop for the plus and negative sense reads
#Going to split the data set into + and - and flip the column names for the negative sense
#This is a bit of an ugly hack, but it should work

positive_sense_recombinations <- subset(recombinations, sense == "+")
nrow(positive_sense_recombinations)
#[1] 4119

negative_sense_recombinations <- subset(recombinations, sense == "-")
nrow(negative_sense_recombinations)
#[1] 4527

#Copy the values into a temporary vector
original_negative_start <- negative_sense_recombinations$start
original_negative_stop <- negative_sense_recombinations$stop

#Save a copy for comparison
negative_sense_recombinations_original <- negative_sense_recombinations

#Now lets replace the values
negative_sense_recombinations$start <- original_negative_stop
negative_sense_recombinations$stop <- original_negative_start

#Looks like that worked, now let's rbind these dataframes together
recombinations_scaled_positive_sense <- rbind(positive_sense_recombinations, negative_sense_recombinations)


#Now let's redraw those plots again!
ggparcoord(subset(recombinations_scaled_positive_sense, reference == "PB2|Segment:1"), columns = 2:3, scale = "globalminmax")
ggparcoord(subset(recombinations_scaled_positive_sense, reference == "PB1|Segment:2"), columns = 2:3, scale = "globalminmax")

#Two adjustments to match Alnaji et al format 
ggparcoord(subset(recombinations_scaled_positive_sense, size > 190 & reference == "PB1|Segment:2" & sample == "barcode95"), columns = c(3, 2), scale = "globalminmax") + coord_flip()
ggparcoord(subset(recombinations_scaled_positive_sense, size > 190 & reference == "PA|Segment:3" & sample == "barcode95"), columns = c(3, 2), scale = "globalminmax") + coord_flip()

#Let's take a peek at the recombination plot now
ggplot(subset(recombinations_scaled_positive_sense, size > 190 & reference == "PB2|Segment:1"), aes(x = stop, y = start, color = number_of_samples)) +
  geom_point(aes(size = count),alpha = 0.6) + 
  scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ treatment)

#Now look at the original plot, but with the rescaled values
ggplot(recombinations_scaled_positive_sense, aes(x = stop, y = start, color = number_of_samples)) +
  geom_point(aes(size = count),alpha = 0.6) + 
  scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference, scales = "free")

#Now look at the original plot, but with the rescaled values, BUT only cannonical
ggplot(subset(recombinations_scaled_positive_sense, size > 190), aes(x = stop, y = start, color = number_of_samples)) +
  geom_point(aes(size = count),alpha = 0.6) + 
  scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference, scales = "free")


#The small deletions are gnawing at me. Can't possibly be so many right?

#Let's look at the size distribution of these small deletions in a actually look at the histogram
ggplot(recombinations, aes(x = size)) + geom_histogram()

#Do the same in log scale
ggplot(recombinations, aes(x = size)) + geom_histogram() + scale_y_log10() + scale_x_log10() + geom_vline(xintercept = c(150, 170, 190), color = "red")

#Cannonicals are definitely different according to sample!!!
ggplot(subset(recombinations, size > 190), aes(x = size)) + geom_histogram() + geom_vline(xintercept = c(150, 170, 190), color = "red") + facet_wrap(~ sample_name)

#Last thing is what is the distribution of only the small deletions
ggplot(subset(recombinations, size < 190), aes(x = size)) + geom_histogram() + geom_vline(xintercept = c(150, 170, 190), color = "red") + facet_wrap(~ sample_name)
#Very interesting, the smaller deletions are way, way more common and they drop off in size
#The pattern is almost identical in each of those samples



##### Are the small deletions real? #####
#After talking with Chris there are a few things to check into to see if the small deletions are real
#1. Check within a given UMI if the 
#2. If the deletions are real, what would be the necessary particle:PFU ratio
#3. Do the short deletions follow the oscillatory pattern
#4. Independently measure the particle:PFU ratio and see if it tracks.


#First the easiest I can do here in R 
#2. Let's caclulate the DVG/SVG ratio.

#Need the files which have read counts
file_names <- NULL
file_names <- list.files("/Users/mixtup/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/output/read_counts/", "*_per_segment_read_count_deduplicated_umi_annotated.txt", full.names = T)


file_names <- file_names[49:96]


per_segment_reads <- NULL
for (file in file_names) {
  print(file)
  df <- read.table(file, sep="\t", quote="")
  df <- cbind(rep(file, nrow(df)),df)
  per_segment_reads <- rbind(per_segment_reads, df)
}

#Name columns
colnames(per_segment_reads) <- c("sample", "reference", "primary_mapped_reads")

#Clean up the filename, remove path
per_segment_reads$sample <- gsub("/Users/mixtup/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/output/read_counts//","", per_segment_reads$sample)
per_segment_reads$sample <- gsub("_per_segment_read_count_deduplicated_umi_annotated.txt","", per_segment_reads$sample)

#Remove the column with no information 
per_segment_reads <- subset(per_segment_reads, reference != "*")

#Adjust the sample column for compatibility in the join
#per_segment$sample <- gsub("bwa_","", per_segment$sample)

#Now summarize from recombination file
#Calculate the number of molecules (corresponding to read becasue) that have been deleted, per segment,  not subsetting by size
per_segment <- recombinations %>% 
  group_by(sample, reference, start, stop, count, size) %>% 
  summarise(unique_reads = n()*count) %>% group_by(sample, reference) %>%  summarise(total_dvg_reads = sum(unique_reads))

#Now join
per_segment <- left_join(per_segment, per_segment_reads, by = c("sample", "reference"))

#And now calculate a DVG to SVG ratio NOTE!!! This calculation has changed as of Oct 11 because
 #I now have only primary reads mapping in the read counts, which means total is the primary_mapped_reads + total_dvg_reads
per_segment <- per_segment %>% 
  mutate(total_reads = primary_mapped_reads+total_dvg_reads) %>%
  mutate(dvg_svg_prop = total_dvg_reads / (primary_mapped_reads+total_dvg_reads)) %>%
  mutate(dvg_svg_perc = dvg_svg_prop*100)

#Now join with data tables
#sample_data <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/minion_genome_sequencing/minion_flu_seq/data/ile_double_umi_runA/sample_data .csv")

#Adjust columns to match for join
colnames(per_segment) <- c("sample", "reference", "total_dvg_reads", "primary_mapped_reads", "total_reads", "dvg_svg_prop", "dvg_svg_perc")

#Join sample
per_segment <- left_join(per_segment, sample_data, by = c("sample"))

#Plot ratios
ggplot(per_segment, aes(x = treatment, y = dvg_svg_perc, color = treatment)) +
  geom_col(position = "dodge") +
  ylab("Percent of DVG Reads") +
  xlab("Segment") +
  theme_classic(base_size = 18) +
  facet_wrap(~reference, scales = "free")

#Get average ratio
per_segment <- per_segment %>% group_by(treatment, sample_name, reference) %>% 
    mutate(average_perc = mean(dvg_svg_perc)) %>%
    mutate(sd_perc = sd(dvg_svg_perc))
  

#Now by average of bioreplicates
ggplot(per_segment, aes(x = treatment, y = average_perc, color = treatment)) +
  geom_col(position = "dodge") +
  ylab("Percent of DVG Reads") +
  xlab("Segment") +
  theme_classic(base_size = 18) +
  facet_wrap(~ reference, scales = "free")

#Flip by Treatment
ggplot(per_segment, aes(x = reference, y = average_perc, color = treatment)) +
  geom_col(position = "dodge") +
  ylab("Percent of DVG Reads") +
  xlab("Segment") +
  theme_classic(base_size = 18) +
  facet_wrap(~ treatment)

#Let's do some sanity checks on DVG and read counts

#What do raw values look like for DVGs 
ggplot(per_segment, aes(x = sample, y = total_dvg_reads, color = treatment)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ reference)
#barcode91 looks low in all of those, but in general all of these are looking 

#What about percentage? 
ggplot(per_segment, aes(x = sample, y = dvg_svg_perc, color = treatment)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ reference)
#lol barcode91 now jumps out as one of the larger ones 

#Let's also look at total reads here
ggplot(per_segment, aes(x = sample, y = total_reads, color = treatment)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ reference)

#In the above, we also need to check how do quality control parameters change this
 #Might require using controls to benchmark the criteria for analyses, using the correlation BUT
  #do we have any technical replicates here? 
  #if not, need to use size and RSC's to make some good guesses as to which are valid

#Next steps are to look at specific deletions:
 #What are common deletions across all samples?
 #Which are deletions shared with DMSO sample?
 #Which are well-supported (RSC, present in replicates), unique deletions in specific treatments

##### More Specific to this paper DVG Drug Screen #####

#Let's take a look at the segment DVG plots and color dots to see if we get any readily visible patters
ggplot(subset(recombinations_scaled_positive_sense, size > 190), aes(x = stop, y = start, color = treatment)) +
  geom_point(aes(size = count, shape = strain),alpha = 0.6) + 
  #scale_size_continuous("Read Count",range = c(0.4, 3)) +
  coord_cartesian() +
  ylab("Donor Site") +
  xlab("Acceptor Site") +
  theme_classic(base_size = 18) +
  #scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference, scales = "free")

#A couple of tests
model0 <- lm(size ~ treatment, data = recombinations_scaled_positive_sense)
summary(model0)
summary.aov(model0)

ggplot(recombinations_scaled_positive_sense, aes(x = treatment, y = size)) +
  geom_point(aes(color = count),alpha = 0.6) + 
  coord_cartesian() +
  ylab("Size of Deletion") +
  xlab("Treatment") +
  theme_classic(base_size = 18) +
  scale_color_gradient("# of Isolates", low = "#87f6ff", high = "#f786ff") +
  facet_wrap(~ reference, scales = "free")




##### / More Specific to this paper DVG Drug Screen #####