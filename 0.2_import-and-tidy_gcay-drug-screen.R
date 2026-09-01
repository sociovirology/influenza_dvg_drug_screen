# this is the script to import the dvg text file, 
#   prep it for binding with the SVG text file, 
#     import the SVG text file,
#      prep it for binding with the DVG text file
#       bind the DVG and SVG files
#        create additional required variables
#         and finally, convert variables to the proper data type (factor, numeric etc)
# BEFORE YOU RUN THIS SCRIPT, 
#     1. move DVG and SVG text files that contain no read info to a different folder
#          these will be the mock-infection controls...barcodes 40, 44, 48, 52, 56, 60
#          parts of the script will error out if it hits an empty file
#          if a non mock-infection Tx group in the DVG file is empty, then there was a ViReMa error
#          if a non mock-infection Tx group in the SVG file is empty, then there was a sequence aligner error
#     2. ensure the import file paths are accurate




library(tidyverse)

#____________________________DVG text file import___________________________


# Get the list of names of the txt files
#file_names_dvg <- list.files("/home/ile/Desktop/play-data/virema_output_drug-screen/dvg_txtfiles", "*_dvg.txt", full.names = T)
file_names_dvg <- list.files("/Users/sociovirology/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/dvg_txtfiles", "*_dvg.txt", full.names = T)

# Import
dvg_txtfile <- NULL
for (file in file_names_dvg) {
  print(file)
  df <- read.table(file, skip=1, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE)
  df <- cbind(df, rep(file, nrow(df)))
  dvg_txtfile <- rbind(dvg_txtfile, df)
}

# Assign variables meaningful names
colnames(dvg_txtfile) <- c("reference" , "reference.acceptor" , "read.name" , "del1" , "del2" , "del3" , "del4" , "del5" , "del6" , "ont.bcd")

# drop observations that contain _NA_ values
dvg_txtfile <- dvg_txtfile %>%
  filter(!is.na(del1))

# create a new variable called "geno.class"
#     if "reference" and "reference.acceptor" are the same, then it's a DVG i.e. same segment recombination
#     if not, then it's a template-switched genome (TSG)
dvg_txtfile <- dvg_txtfile %>%
  mutate(geno.class = ifelse(reference == reference.acceptor, "DVG", "TSG")
  )

#_________________filter-out un-collapsed recombination events_______________________
# the dataset contains collapsed DVG records, but also uncollapsed recombination events
#        the template-switched and hybrids (deletion + template-switched) chimeric recombinations that share the same read name
#        so FOR NOW, we filter them out
dvg_txtfile <- dvg_txtfile %>%
  filter(!(duplicated(read.name) | duplicated(read.name, fromLast = TRUE)))

#___________________filter out DVGs that contain more than 2x deletions__________________
# it's a current glitch in ViReMa I'm still TSing
# the sum of their deletions is greater than their genomic length, which is clearly wrong
dvg_txtfile <- dvg_txtfile %>%
  filter( del3 == 0 & del4 == 0 & del5 == 0 & del6 == 0)

# Extract the sample name i.e. ONT barcode# from the filepath variable
# We want the "sample" aka "ont.bcd"  variable to say what barcode the given dvg event stems frm
# The "ont.bcd" variable values have that info, but only as part of a super bulky filepath
# The nxt two cmds successively trim the filepath until only the	_barcode##_ is rmg
#dvg_txtfile$ont.bcd <- gsub("/home/ile/Desktop/play-data/virema_output_drug-screen/dvg_txtfiles/","", dvg_txtfile$ont.bcd)
dvg_txtfile$ont.bcd <- gsub("/Users/sociovirology/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/dvg_txtfiles/","", dvg_txtfile$ont.bcd)
dvg_txtfile$ont.bcd <- gsub("_dvg.txt","", dvg_txtfile$ont.bcd)


#
#
# prepping dvg_txtfile df for the binding with svg_txtfile
#     this step will output a new object called "dvg.bind"
#
# the current variables names in dvg_txtfile are:
# "reference" , "reference.acceptor" , "read.name" , "del1" , "del2" , "del3" , "del4" , "del5" , "del6" , "ont.bcd")
#
# The only starter variables I need from the	dvg_txtfiles	dataframe are:
#	"reference" , "reference.acceptor" , "read.name" , "del1" , "del2", "geno.class" , and "ont.bcd"
#
# I'm also adding	filter(!is.na(del1))	in order to collect only records w/o _NA_ ...just a safety redundancy...so when I fuse dvg.bind with svg.bind, there will be no _NA_ cells , only _0_ 

dvg.bind <- dvg_txtfile %>%
  dplyr::select(reference, reference.acceptor , read.name , del1 , del2, geno.class , ont.bcd)  %>%
  filter(!is.na(del1))





#
#
#____________________________________SVG text file import_________________________________________

# Get the list of names of the SVG txt files
file_names_svg <- list.files("/Users/sociovirology/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/svg_txtfiles", "*_svg.txt", full.names = T)

# Import
svg_txtfile <- NULL
for (file in file_names_svg) {
  print(file)
  df <- read.table(file, skip=1, sep="\t", stringsAsFactors=FALSE, quote="")
  df <- cbind(df, rep(file, nrow(df)))
  svg_txtfile <- rbind(svg_txtfile, df)
}

# Rename columns more meaningfully
colnames(svg_txtfile) <- c("read.name", "reference", "read.lgt", "geno.class", "ont.bcd")

# Extract the sample name i.e. ONT barcode# from the filepath variable
# We want the "sample" aka "ont.bcd"  variable to say what barcode the given dvg event stems frm
# The "ont.bcd" variable values have that info, but only as part of a super bulky filepath
# The nxt two cmds successively trim the filepath until only the	_barcode##_ is rmg
svg_txtfile$ont.bcd <- gsub("/Users/sociovirology/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/svg_txtfiles/","", svg_txtfile$ont.bcd)
svg_txtfile$ont.bcd <- gsub("_svg.txt","", svg_txtfile$ont.bcd)


#
#
# prepping svg_txtfile df for the binding with dvg_txtfile
#     this step will output a new object called "svg.bind"
#
# The both dataframes currently have the following variables:

# dvg.bind:		"reference" , "reference.acceptor" , "read.name" , "del1" , "del2", "geno.class" , and "ont.bcd"
# svg_txtfile:		"read.name", "reference", "read.lgt", "geno.class", and "ont.bcd"


# svg.bind transformation will be as follows:
#	starting w/	svg_txtfile	
#	add the following variables:	"reference.acceptor", "del1" , and "del2"
#	then select() the same ones as are in dvg.bind:	"reference" , "reference.acceptor" , "read.name" , "del1" , "del2", "geno.class"


svg.bind <- svg_txtfile %>%
  mutate(del1 = 0 ) %>%
  mutate(del2 = 0 ) %>%
  mutate(reference.acceptor = NA) %>%
  dplyr::select(reference, reference.acceptor , read.name , del1 , del2, geno.class , ont.bcd)




#____________________DVG-SVG bind_______________________________
# Now, we bind the DVG and SVG dataframes into a single object so we can calculate:
#     i. total viral genomes (TVG)  
#     ii. relative abundance or proportion of DVGs (dvg/tvg)... "prop.dvg" variable

# This is where I bind the DVG and SVG dataframes to gen "count.dvg" and "prop.DVG" variables

# GIVEN:
#	a dataframe called	dvg.bind	w/the ff variables:	"reference" , "reference.acceptor" , "read.name" , "del1" , "del2", "geno.class" , and "ont.bcd"
#	a dataframe called	svg.bind	w/the ff variables:	"reference" , "reference.acceptor" , "read.name" , "del1" , "del2", "geno.class" , and "ont.bcd"


# TASK:
#	Bind them, then create the new variables necessary for downstream analysis
#
# bind dvg.bind and svg.bind atop each other
dvg.svg.bind <- bind_rows(dvg.bind, svg.bind)

# Before I create new variables, I use str() to confirm the source variable are the right data-type
#
# > str(dvg.svg.bind)
# 'data.frame':	4775515 obs. of  7 variables:
# $ reference         : chr  "PB1|Segment:2_RevStrand" "PB2|Segment:1_RevStrand" "M|Segment:7" "PA|Segment:3" ...
# $ reference.acceptor: chr  "PB1|Segment:2_RevStrand" "PB2|Segment:1_RevStrand" "M|Segment:7" "PA|Segment:3" ...
# $ read.name         : chr  "16faf580-5ec3-5a17-a9f4-42d803852339_AATTGGATAGTG_Fuzz=" "d915db7f-9c31-4df4-8809-df3660aeb3dd_TTTTCGTCTGTC_Fuzz=" "5799f83e-8b27-4769-8d25-69250ca86167_TAGTAAATGTTT_Fuzz=" "afdc1433-74ee-42dc-960e-48154f77dfcb_GTCGGTCCGTAA_Fuzz=" ...
# $ del1              : num  1635 1514 142 1350 229 ...
# $ del2              : num  0 0 0 0 0 0 0 0 0 0 ...
# $ geno.class        : chr  "DVG" "DVG" "DVG" "DVG" ...
# $ ont.bcd           : chr  "barcode01" "barcode01" "barcode01" "barcode01" ...
#
#
# Looks like I need to
#	make reference a factor variable
#	make reference.acceptor a factor variable
#	make geno.class a factor variable
#	make ont.bcd a factor variable

dvg.svg.bind <- dvg.svg.bind %>% 
  mutate(reference=as.factor(reference)) %>% 
  mutate(reference.acceptor=as.factor(reference.acceptor)) %>% 
  mutate(ont.bcd=as.factor(ont.bcd)) %>% 
  mutate(geno.class=as.factor(geno.class))


# Now I derive new variables from existing ones, and re-level the new variables as necessary
# Here's what I need to do
#	create "geno.sgmt" from data in "reference" i.e. PB2|Segment:1 in2 PB2	...then relevel it	PB2 > PB1 > PA > HA > NP > NA > M > NS
#	create "geno.sgmt.acceptor" from data in "reference.acceptor"	         	...then relevel it	PB2 > PB1 > PA > HA > NP > NA > M > NS
#	create "bioreplicate" from data in "ont.bcd" 			                    	...then relevel it	1 > 2 > 3
#	create "Tx" from data in "ont.bcd" 				                            	...then relevel it	dH2O > DMSO > Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
#	create "veh.Tx" from data in "Tx"     		                	         		...then relevel it	DMSO > dH2O
# create "grouping.Tx" from data in "Tx"                                     ...then relevel it  DMSO > dH2O > prelim
#	create "host.ptg" from data in "ont.bcd"			                        	...then relevel it	MDCK.mock > MDCK.CA09(H1N1) > MDCK.TX12(H3N2)


#_________________create the new variable "geno.sgmt" from data in "reference" variable
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(geno.sgmt = case_when(reference %in% c("PB2|Segment:1" , "PB2|Segment:1_RevStrand") ~ "PB2",
                               reference %in% c("PB1|Segment:2" , "PB1|Segment:2_RevStrand") ~ "PB1",
                               reference %in% c("PA|Segment:3" , "PA|Segment:3_RevStrand") ~ "PA",
                               reference %in% c("HA|Segment:4" , "HA|Segment:4_RevStrand") ~ "HA",
                               reference %in% c("NP|Segment:5" , "NP|Segment:5_RevStrand") ~ "NP",
                               reference %in% c("NA|Segment:6" , "NA|Segment:6_RevStrand") ~ "NA",
                               reference %in% c("M|Segment:7" , "M|Segment:7_RevStrand") ~ "M",
                               reference %in% c("NS|Segment:8" , "NS|Segment:8_RevStrand") ~ "NS" )
  )


#_________________create the new variable "geno.sgmt.acceptor" from data in "reference.acceptor" variable
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(geno.sgmt.acceptor = case_when(reference.acceptor %in% c("PB2|Segment:1" , "PB2|Segment:1_RevStrand") ~ "PB2",
                                        reference.acceptor %in% c("PB1|Segment:2" , "PB1|Segment:2_RevStrand") ~ "PB1",
                                        reference.acceptor %in% c("PA|Segment:3" , "PA|Segment:3_RevStrand") ~ "PA",
                                        reference.acceptor %in% c("HA|Segment:4" , "HA|Segment:4_RevStrand") ~ "HA",
                                        reference.acceptor %in% c("NP|Segment:5" , "NP|Segment:5_RevStrand") ~ "NP",
                                        reference.acceptor %in% c("NA|Segment:6" , "NA|Segment:6_RevStrand") ~ "NA",
                                        reference.acceptor %in% c("M|Segment:7" , "M|Segment:7_RevStrand") ~ "M",
                                        reference.acceptor %in% c("NS|Segment:8" , "NS|Segment:8_RevStrand") ~ "NS" )
  )


#__________________create the new variable "bioreplicate" from data in "ont.bcd" variable
#	bioreplicates are the #s 1 , 2 , 3
# the following barcodes belong to their assigned bioreplicates, but as they are mock infection controls they actually have no data on the board...   40 44 48 52 56 60
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(bioreplicate = case_when(
    ont.bcd %in% c("barcode01", "barcode02", "barcode03", "barcode04", "barcode13", "barcode14", "barcode15", "barcode16", "barcode25", "barcode26", "barcode27", "barcode28", "barcode37", "barcode38", "barcode39", "barcode40", "barcode49", "barcode50", "barcode51", "barcode52", "barcode61", "barcode62", "barcode63", "barcode64", "barcode73", "barcode74", "barcode75", "barcode76", "barcode85", "barcode86", "barcode87", "barcode88") ~ "1",
    ont.bcd %in% c("barcode05", "barcode06", "barcode07", "barcode08", "barcode17", "barcode18", "barcode19", "barcode20", "barcode29", "barcode30", "barcode31", "barcode32", "barcode41", "barcode42", "barcode43", "barcode44", "barcode53", "barcode54", "barcode55", "barcode56", "barcode65", "barcode66", "barcode67", "barcode68", "barcode77", "barcode78", "barcode79", "barcode80", "barcode89", "barcode90", "barcode91", "barcode92") ~ "2",
    ont.bcd %in% c("barcode09", "barcode10", "barcode11", "barcode12", "barcode21", "barcode22", "barcode23", "barcode24", "barcode33", "barcode34", "barcode35", "barcode36", "barcode45", "barcode46", "barcode47", "barcode48", "barcode57", "barcode58", "barcode59", "barcode60", "barcode69", "barcode70", "barcode71", "barcode72", "barcode81", "barcode82", "barcode83", "barcode84", "barcode93", "barcode94", "barcode95", "barcode96") ~ "3" )
  ) 


#____________________create the new variable "Tx" from data in "ont.bcd" variable
#	Tx's are dH2O > DMSO > Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(Tx = case_when(
    ont.bcd %in% c("barcode01", "barcode49", "barcode05", "barcode53", "barcode09", "barcode57") ~ "Lepto",
    ont.bcd %in% c("barcode02", "barcode50", "barcode06", "barcode54", "barcode10","barcode58") ~ "Tolyp",
    ont.bcd %in% c("barcode03", "barcode51", "barcode07", "barcode55", "barcode11", "barcode59") ~ "Uri",
    ont.bcd %in% c("barcode04", "barcode64", "barcode08","barcode68", "barcode12", "barcode72") ~ "Favp",
    ont.bcd %in% c("barcode25", "barcode40", "barcode73", "barcode29", "barcode44", "barcode77", "barcode33", "barcode48", "barcode81") ~ "dH2O",
    ont.bcd %in% c("barcode26", "barcode74", "barcode30", "barcode78", "barcode34", "barcode82") ~ "Oscil",
    ont.bcd %in% c("barcode27", "barcode75", "barcode31", "barcode79", "barcode35", "barcode83") ~ "Insu",
    ont.bcd %in% c("barcode28", "barcode88", "barcode32", "barcode92", "barcode36", "barcode96") ~ "4-OI",
    ont.bcd %in% c("barcode37", "barcode52", "barcode85", "barcode41", "barcode56", "barcode89", "barcode45", "barcode60", "barcode93") ~ "DMSO",
    ont.bcd %in% c("barcode38", "barcode86", "barcode42", "barcode90", "barcode46", "barcode94") ~ "Nosto",
    ont.bcd %in% c("barcode39", "barcode87", "barcode43", "barcode91", "barcode47", "barcode95") ~ "MK2206",
    ont.bcd %in% c("barcode16", "barcode76", "barcode20", "barcode80", "barcode24", "barcode84") ~ "UK5099",
    ont.bcd %in% c("barcode13", "barcode61", "barcode17", "barcode65", "barcode21", "barcode69") ~ "Alpe",
    ont.bcd %in% c("barcode14", "barcode62", "barcode18", "barcode66", "barcode22", "barcode70") ~ "Synec",
    ont.bcd %in% c("barcode15", "barcode63", "barcode19", "barcode67", "barcode23", "barcode71") ~ "Ado" )
  )

#_________________create the new variable "veh.Tx" from data in "Tx" variable
#	DMSO vehicle:	"DMSO" , "Alpe" , "MK2206" , "4-OI" , "UK5099" , "Favp"
#	dH2O vehicle:	"dH2O" , "Insu" , "Ado" , "Uri" , "Lepto" , "Nosto" , "Oscil" , "Synec" , "Tolyp"
# this is to allow me easily filter Tx groups by their vehicle during statistics and graphing
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(veh.Tx = case_when(
    Tx %in% c("DMSO", "Alpe", "MK2206", "4-OI", "UK5099", "Favp") ~ "DMSO",
    Tx %in% c("dH2O", "Insu", "Ado", "Uri", "Lepto", "Nosto", "Oscil", "Synec", "Tolyp" ) ~ "dH2O")
  )


#__________________create the new variable "grouping.Tx" from data in "Tx" variable
# the groupings for statistical analysis are "DMSO" , "dH2O" , "prelim"
# this will allow me efficiently filter Tx's that will undergo ANOVA/MANOVA/PERMANOVA together
#     metabolites + metabolic signaling molecules are grouped by their vehicle; DMSO or dH2O
#     whereas, my exploratory Tx's will be compared against e/o w/o a vehicle-ctrl...hence "prelim" for prelim data
#	DMSO vehicle:	"DMSO" , "Alpe" , "MK2206" , "4-OI" , "UK5099"
#	dH2O vehicle:	"dH2O" , "Insu" , "Ado" , "Uri"
# prelim/exploratory Tx: "Lepto" , "Nosto" , "Oscil" , "Synec" , "Tolyp" , "Favp"
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(grouping.Tx = case_when(
    Tx %in% c("DMSO", "Alpe", "MK2206", "4-OI", "UK5099") ~ "DMSO",
    Tx %in% c("dH2O", "Insu", "Ado", "Uri") ~ "dH2O",
    Tx %in% c("Lepto", "Nosto", "Oscil", "Synec", "Tolyp" , "Favp") ~ "prelim")
  )


#___________________create the new variable "host.ptg" from data in "ont.bcd" variable
#	Inoculum choices are MDCK.CA09(H1N1) , MDCK.TX12(H3N2)
# there's no 'MDCK.mock'...those barcodes had no data and need to be removed to prevent an error during -R- import
dvg.svg.bind <- dvg.svg.bind %>%
  mutate(host.ptg = case_when(
    ont.bcd %in% c("barcode01", "barcode02", "barcode03", "barcode04", "barcode05", "barcode06", "barcode07", "barcode08", "barcode09", "barcode10", "barcode11", "barcode12", "barcode13", "barcode14", "barcode15", "barcode16", "barcode17", "barcode18", "barcode19", "barcode20", "barcode21", "barcode22", "barcode23", "barcode24", "barcode25", "barcode26", "barcode27", "barcode28", "barcode29", "barcode30", "barcode31", "barcode32", "barcode33", "barcode34", "barcode35", "barcode36", "barcode37", "barcode38", "barcode39", "barcode41", "barcode42", "barcode43", "barcode45", "barcode46", "barcode47") ~ "MDCK.TX12(H3N2)",
    ont.bcd %in% c("barcode49", "barcode50", "barcode51", "barcode53", "barcode54", "barcode55", "barcode57", "barcode58", "barcode59", "barcode61", "barcode62", "barcode63", "barcode64", "barcode65", "barcode66", "barcode67", "barcode68", "barcode69", "barcode70", "barcode71", "barcode72", "barcode73", "barcode74", "barcode75", "barcode76", "barcode77", "barcode78", "barcode79", "barcode80", "barcode81", "barcode82", "barcode83", "barcode84", "barcode85", "barcode86", "barcode87", "barcode88", "barcode89", "barcode90", "barcode91", "barcode92", "barcode93", "barcode94", "barcode95", "barcode96") ~ "MDCK.CA09(H1N1)" )
  ) 


#___________________Now I relevel all the factor variables accordingly
# use str() to see which variables need data-type adjustments
# str(dvg.svg.bind)
#
# 'data.frame':	4775515 obs. of  14 variables:
# $ reference         : Factor w/ 16 levels "HA|Segment:4",..: 14 16 3 11 10 9 3 12 12 16 ...
# $ reference.acceptor: Factor w/ 16 levels "HA|Segment:4",..: 14 16 3 11 10 11 3 12 12 16 ...
# $ read.name         : chr  "16faf580-5ec3-5a17-a9f4-42d803852339_AATTGGATAGTG_Fuzz=" "d915db7f-9c31-4df4-8809-df3660aeb3dd_TTTTCGTCTGTC_Fuzz=" "5799f83e-8b27-4769-8d25-69250ca86167_TAGTAAATGTTT_Fuzz=" "afdc1433-74ee-42dc-960e-48154f77dfcb_GTCGGTCCGTAA_Fuzz=" ...
# $ del1              : num  1635 1514 142 1350 229 ...
# $ del2              : num  0 0 0 0 0 0 0 0 0 0 ...
# $ geno.class        : Factor w/ 3 levels "DVG","SVG","TSG": 1 1 1 1 1 3 1 1 1 1 ...
# $ ont.bcd           : Factor w/ 90 levels "barcode01","barcode02",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ geno.sgmt         : chr  "PB1" "PB2" "M" "PA" ...
# $ geno.sgmt.acceptor: chr  "PB1" "PB2" "M" "PA" ...
# $ bioreplicate      : chr  "1" "1" "1" "1" ...
# $ Tx                : chr  "Lepto" "Lepto" "Lepto" "Lepto" ...
# $ veh.Tx            : chr  "dH2O" "dH2O" "dH2O" "dH2O" ...
# $ grouping.Tx       : chr  "prelim" "prelim" "prelim" "prelim" ...
# $ host.ptg          : chr  "MDCK.TX12(H3N2)" "MDCK.TX12(H3N2)" "MDCK.TX12(H3N2)" "MDCK.TX12(H3N2)" ...

# looks like I need to as.factor() the new variables first
dvg.svg.bind <- dvg.svg.bind %>% 
  mutate(geno.sgmt=as.factor(geno.sgmt)) %>% 
  mutate(geno.sgmt.acceptor=as.factor(geno.sgmt.acceptor)) %>% 
  mutate(bioreplicate=as.factor(bioreplicate)) %>% 
  mutate(Tx=as.factor(Tx)) %>% 
  mutate(veh.Tx=as.factor(veh.Tx)) %>% 
  mutate(grouping.Tx=as.factor(grouping.Tx)) %>% 
  mutate(host.ptg=as.factor(host.ptg))


# Now to relevel...str() shows that below variables need re-leveling as follows:
#	"geno.sgmt"	             	PB2 > PB1 > PA > HA > NP > NA > M > NS
#	"geno.sgmt.acceptor"     	PB2 > PB1 > PA > HA > NP > NA > M > NS
#	"Tx"		                 	dH2O > DMSO > Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
#	"veh.Tx"	              	DMSO > dH2O
# "grouping.Tx"                DMSO > dH2O > prelim
#	"host.ptg"		MDCK.mock > MDCK.CA09(H1N1) > MDCK.TX12(H3N2)

# library(forcats) ...tidyverse already loaded
dvg.svg.bind$geno.sgmt <- fct_relevel(dvg.svg.bind$geno.sgmt, "PB2" , "PB1" , "PA" , "HA" , "NP" , "NA" , "M" , "NS")
dvg.svg.bind$geno.sgmt.acceptor <- fct_relevel(dvg.svg.bind$geno.sgmt.acceptor, "PB2" , "PB1" , "PA" , "HA" , "NP" , "NA" , "M" , "NS")
dvg.svg.bind$Tx <- fct_relevel(dvg.svg.bind$Tx, "dH2O", "DMSO", "Insu", "Alpe", "MK2206", "4-OI", "UK5099", "Ado", "Uri", "Favp", "Lepto", "Nosto", "Oscil", "Synec", "Tolyp")
dvg.svg.bind$veh.Tx <- fct_relevel(dvg.svg.bind$veh.Tx, "DMSO", "dH2O")
dvg.svg.bind$grouping.Tx <- fct_relevel(dvg.svg.bind$grouping.Tx, "DMSO", "dH2O", "prelim")
dvg.svg.bind$host.ptg <- fct_relevel(dvg.svg.bind$host.ptg, "MDCK.CA09(H1N1)" , "MDCK.TX12(H3N2)")



#test with ... levels(dvg.svg.bind$Tx)
# use str(dvg.svg.bind) to ensure the dvg.svg.bind dataframe is fully set up

# free up room in R environment
rm("df" , "dvg_txtfile" , "dvg.bind", "file" , "file_names_dvg" , "file_names_svg" , "svg_txtfile", "svg.bind")


#_____________dvg.svg.bind explained:

# each observation is a single Influenza A genomic segment...
#     DVG or SVG, 
#     from parent flu stock,  
#     or the progeny stock from the alpelisib dosing experiment

# the variables
# reference         : the sequence alignment identity of the read...for DVGs, this is the recombination donor segment
# reference.acceptor: the sequence alignment identity of the DVG recombination acceptor...'NA' for SVGs
# read.name         : the read name from the fastq
# del1              : the first deletion-type recombination event observed in the genomic segment; reading the gene sequence from left to right (0 to n)
# del2              : the second
# geno.class        : is the read a DVG or SVG
# ont.bcd           : the Oxford Nanopore (ONT) barcode#
# geno.sgmt         : same as 'reference', but in a coding friendly format
# geno.sgmt.acceptor: same as 'reference.acceptor', but in a coding friendly format
# bioreplicate      : bioreplicates of the drug-dosing screen...three experiments were run on three different days
# Tx                : Each administered drug, metabolite, or biologic
# veh.Tx            : the treatments grouped by their vehicle...DMSO or dH2O
# grouping.Tx       : the treatments grouped by their similarity...
# host.ptg          : host cell, and the pathogen it was inoculated with
#
#
##	DMSO vehicle:	"DMSO" , "Alpe" , "MK2206" , "4-OI" , "UK5099"
#	dH2O vehicle:	"dH2O" , "Insu" , "Ado" , "Uri"
# prelim/exploratory Tx: "Lepto" , "Nosto" , "Oscil" , "Synec" , "Tolyp" , "Favp"

#
#
#