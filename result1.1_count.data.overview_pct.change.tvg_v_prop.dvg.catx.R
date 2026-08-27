# Drug screen genome counting assay (gcay) results
#
# _______Results_____________
# RESULT 1.2: is the plot of 
#     (i)total viral genomes and 
#     (ii)proportion of defective viral genomes 
# for each Tx
#
# but first, let's look at the data



library(tidyverse)

#____FIRST: 
# summarise total genomes (count.tvg), 
#           total dvgs (count.dvg), and 
#           relative abundance--aka proportion--of DVGs (prop.dvg)
gcay.screen.count <- dvg.svg.bind %>%
  group_by(host.ptg, Tx, veh.Tx, bioreplicate) %>%
  summarise(count.tvg = n(),
            count.dvg = sum(geno.class == "DVG"),
            prop.dvg = (sum(geno.class == "DVG")) / n()    ) %>%
  ungroup() 


#___SECOND:
# create an object of averaged bioreplicates
gcay_screen_bio.avg <-
  gcay.screen.count %>%
  group_by(host.ptg, Tx, veh.Tx) %>%
  summarise(n = n(),
            mean.count.tvg = mean(count.tvg),
            sd.tvg = sd(count.tvg),
            sem.tvg = sd.tvg /sqrt(n()),
            
            mean.count.dvg = mean(count.dvg),
            sd.dvg = sd(count.dvg),
            sem.dvg = sd.dvg /sqrt(n()),
            
            mean.prop.dvg = mean(prop.dvg),
            sd.prop.dvg = sd(prop.dvg),
            sem.prop.dvg = sd.prop.dvg /sqrt(n())
  ) %>%
  ungroup() 


#_________Now we can overview the data________


# Visualize the full spread of total genomes produced across the bioreplicates
overview.1.gcay.data.spread <- 
  gcay.screen.count %>%
  ggplot(aes(x=Tx, y=count.tvg, color=bioreplicate)) +
  geom_jitter(stat = "identity", width = 0.1, height = 0 ) +
  labs(x = "Treatment", y = "TVG Count") +
  #geom_text(aes(label = count.tvg) , size = 2.4 , vjust = -0.5 , hjust = 1.15) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "bottom" ) +
  facet_wrap(~host.ptg)

# data spread is high and low across the strains and treatments
#    best to normalize by deriving the percent.change from each bioreplicate's vehicle treatment groups (DMSO + dH2O)


#____Visualize the averaged total genomes from the bioreplicates (n=3)
overview.2.gcay.data.spread <-
  gcay_screen_bio.avg %>%
  ggplot( aes(x=Tx, y=mean.count.tvg, fill=host.ptg)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatment", y = "Mean Count of TVGs") +
  geom_text(aes(label = round(mean.count.tvg, digits = 0) ) , color="black", 
            size = 2 , angle=45 , vjust = -2 , hjust = 0) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "none" ) +
  geom_errorbar(aes(x=Tx,               
                    ymin=mean.count.tvg - sem.tvg , 
                    ymax=mean.count.tvg + sem.tvg , 
                    width=0.1)
                ) +
  facet_wrap(~host.ptg)

# barchart and error bars do a much better job of visually communicating the variance
#     as already mentioned, normalizing by percent.change is the way to go


#____Visualize the averaged proportion of DVGs spawned under each Tx
overview.4.gcay.data.spread <-
  gcay_screen_bio.avg %>%
  ggplot( aes(x=Tx, y=mean.prop.dvg, fill=host.ptg)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatment", y = "Mean Proportion of DVGs") +
  geom_text(aes(label = round(mean.prop.dvg, digits = 3) ) , color="black", 
            size = 2 , angle=45 , vjust = -2 , hjust = 0) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "none" ) +
  geom_errorbar(aes(x=Tx,               
                    ymin=mean.prop.dvg - sem.prop.dvg , 
                    ymax=mean.prop.dvg + sem.prop.dvg , 
                    width=0.1)
  ) +
  facet_wrap(~host.ptg)

# We expected the proportion of deletion-type recombination genomes (DVGs) to still be very low at 18 h.p.i. ... 
#      this was the case for CA09, with a range between 0.019-0.024
#      but not for TX12, which had a more pronounced range of 0.023-0.035
# As with the Alpelisib study, we can overcome this 'limit of detection' problem in two ways:
#     i.  Using %change will allow us better capture the individuality of the bioreplicates,
#     ii. going from the 'total genome' level down to the 'per segment' level ... there is likely variation in each segment's response to treatment ... they behave as individuals in some regards e.g. reassortment
#
# FINALLY, using %change for prop.dvg and total genomes will allow us put both variables on the same plot in a way that makes comparison intuitive
# data spread is low in most cases, which is good.




#_______export the plots...conveniently export plots w/ggsave()
# Set desired save directory
save_directory <- "~/R_gcay_screen"
# save to above directory
ggsave(file.path(save_directory, "count.tvg.data.spread.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.1.gcay.data.spread, 
       width = 7 , height = 4)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "mean.count.tvg.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.2.gcay.data.spread, 
       width = 7 , height = 3)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "mean.prop.dvg.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.4.gcay.data.spread, 
       width = 7 , height = 3)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements


# In light of 
#    i. large spread of the data at the level of total genomes, and
#    ii. the limit of detection for change to prop.dvg, also at the level of total genomes
#     it is more sound to conduct analysis at the per segment level



#______THIRD: 
# summarise per-segment count of total genomes (count.tvg.sgmt), 
#                                total dvgs (count.dvg.sgmt), and 
#                                relative abundance--aka proportion--of DVGs (prop.dvg.sgmt)
gcay.screen.count.sgmt <- 
  dvg.svg.bind %>%
  group_by(host.ptg, Tx, veh.Tx, bioreplicate, geno.sgmt) %>%
  summarise(count.tvg.sgmt = n(),
            count.dvg.sgmt = sum(geno.class == "DVG"),
            prop.dvg.sgmt = (sum(geno.class == "DVG")) / n()    ) %>%
  ungroup()


#______FOURTH:
# create an object of averaged bioreplicates from the above object of averaged tech replicates
gcay_screen_bio.avg.sgmt <-
  gcay.screen.count.sgmt %>%
  group_by(host.ptg, Tx, veh.Tx, geno.sgmt) %>%
  summarise(n = n(),
            mean.count.tvg.sgmt = mean(count.tvg.sgmt),
            sd.tvg.sgmt = sd(count.tvg.sgmt),
            sem.tvg.sgmt = sd.tvg.sgmt /sqrt(n()),
            
            mean.count.dvg.sgmt = mean(count.dvg.sgmt),
            sd.dvg.sgmt = sd(count.dvg.sgmt),
            sem.dvg.sgmt = sd.dvg.sgmt /sqrt(n()),
            
            mean.prop.dvg.sgmt = mean(prop.dvg.sgmt),
            sd.prop.dvg.sgmt = sd(prop.dvg.sgmt),
            sem.prop.dvg.sgmt = sd.prop.dvg.sgmt /sqrt(n())
  ) %>%
  ungroup() 



#_________Now we can overview the per segment data________


#____Visualize the full spread of the per-segment break-out of total genomes produced across the bioreplicates
overview.1.gcay.data.spread.sgmt <- 
  gcay.screen.count.sgmt %>%
  ggplot(aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_jitter(stat = "identity", width = 0.1, height = 0, alpha = 0.7) +
  labs(x = "Treatment", y = "TVG Count per Segment") +
  #geom_text(aes(label = count.tvg.sgmt) , size = 2.4 , vjust = -0.5 , hjust = 1.15) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "bottom" ) +
  facet_grid(geno.sgmt~host.ptg,
             scales = "free_y")
# very promising. As with the Alpelisib study, there is a clearer sense of unique responses at the segment level 
#   let's take the mean and standard error of total genomes to better visualize variance


#___Visualize the per segment mean total genomes (n=3 bioreplicates)
overview.2.gcay.data.spread.sgmt <- 
  gcay_screen_bio.avg.sgmt %>%
  ggplot(aes(x=Tx, y=mean.count.tvg.sgmt, fill=host.ptg)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatment", y = "Mean TVG Count per Segment") +
  geom_text(aes(label = round(mean.count.tvg.sgmt, digits = 0) ) , 
            color="black", size = 2 , angle=45 , vjust = -2 , hjust = 0) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "none" ) +
  geom_errorbar(aes(x=Tx,              
                    ymin=mean.count.tvg.sgmt - sem.tvg.sgmt , 
                    ymax=mean.count.tvg.sgmt + sem.tvg.sgmt , 
                    width=0.1)
  ) +
  facet_grid(geno.sgmt~host.ptg,
             #scales = "free_y"           # ALERT: Don't "free_y" on this plot; the look is quite misleading. The shared axis is more informative
             )
# The unique response of the test groups relative to vehicle is more evident at the segment level
#     Leptolyngya in particular appears to decimate population of HA and NA segments of TX12, while causing every other segment to increase in total counts
# the per-segment level may be the minimum resolution from where to do analysis



#____Visualize the mean proportion of each segment’s total genomes that are DVGs (n=3 bioreplicates)
#        i.e. each segment's DVG-spawning probability in the form of a fraction
#              e.g. the proportion of total PB2 genomes that are DVGs, and so on for the remaining segments
overview.4.gcay.data.spread.sgmt <- 
  gcay_screen_bio.avg.sgmt %>%
  ggplot(aes(x=Tx, y=mean.prop.dvg.sgmt, fill=host.ptg)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatment", y = "Mean DVG Proportion per Segment") +
  geom_text(aes(label = round(mean.prop.dvg.sgmt, digits = 3) ) , 
            color="black", size = 2 , angle=45 , vjust = -2 , hjust = 0) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "none" ) +
  geom_errorbar(aes(x=Tx,              
                    ymin=mean.prop.dvg.sgmt - sem.prop.dvg.sgmt , 
                    ymax=mean.prop.dvg.sgmt + sem.prop.dvg.sgmt , 
                    width=0.1)
  ) +
  facet_grid(geno.sgmt~host.ptg,
             #scales = "free_y"           # ALERT: Don't "free_y" on this plot; the look is quite misleading. The shared axis is more informative
  )
# a very telling plot. A few things are clear about strain-specific outcomes in DVG proportion
#
# REMINDER: Proportions here means the proportion of a segment's total genomes that are DVGs
#   So you can eyeball the bars from   (i)  left to right..."how does each Tx affect each segment's probability of generating a DVG at 18 h.p.i. ?" ... It's clear that Insulin and Adenosine have a strong positive effect on the DVG-spawning probability of the IAVpol segments for both CA09 and TX12
#                            or from   (ii) top to bottom..."in a given Tx at 18hpi, what is each segment's probability of spawning a DVG ?" ... In all Tx's, for both CA09 and TX12, the IAVpol segments appear to have the highest DVG-spawning probability ... Outcomes become more randomized for the remaining TX12 segments, but for CA09, the HA--and then NP--segments have the next largest DVG-spawning probability ... it is from here that outcomes become more randomized for the remaining TX12 segments
#
#
# these findings beg the question of the total population/community of DVGs, and the per segment breakdown




#____Visualize the proportion of each segment's DVGs as a % of total dvg's (n=3 bioreplicates)
#   i.e. of all DVG's in a sample, what proportion is PB2, PB1, and so on for the remaining segments
#    said another way, for PB2, if you randomly pull a DVG from that sample, what is the probability that it will be a PB2 DVG?
#
# use a grouped-%stacked barchart where:
# Group is 'Tx' 'Alpe.uM.fctr'
# Sub-group is 'bioreplicate'
# Sub-sub-group is 'geno.sgmt'
overview.5.gcay.data.spread.sgmt <-
  dvg.svg.bind %>%
  filter(host.ptg %in% c("MDCK.CA09(H1N1)", "MDCK.TX12(H3N2)"),
         geno.class=="DVG") %>%
  group_by(host.ptg, Tx, bioreplicate, geno.sgmt) %>%
  summarise(count.dvg = n()  ) %>%
  mutate(pct.dvg = signif( count.dvg / sum(count.dvg) *100 , digits = 3)
  ) %>%
  ungroup() %>%                     # not quite sure why, but you need this ungroup() for follow-up factor re-leveling to work
  mutate(interaction_group = factor(interaction(Tx, bioreplicate) )
  ) %>%
  mutate(interaction_group = fct_inorder(interaction_group)
  ) %>%
  ggplot(aes(x = interaction_group, y = pct.dvg, fill = geno.sgmt)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Treatments by bioreplicate", y = "Relative Abundance (as percentage of Total DVGs)") +
  geom_text(aes(label = count.dvg) , size = 2.4, angle=69  , position = position_stack(vjust = 0.5)) +  # this prints the count.dvg value onto each geno.sgmt
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7)) +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = margin(c(0, 0, 8, 0)),    # Adjust the legend margins (top, right, bottom, left)
        legend.key.height = unit(0.1, "cm")        # Adjust the height of legend keys to control the number of rows
  ) +
  guides(fill = guide_legend(ncol = 8)) +          # Set the number of columns in the legend...need this to flatten legend from 2x4 grid to 1x8 grid
  labs(fill = "genomic segment") +                 # insert your custom legend title
  theme(plot.margin = margin(t = -1, r = 10, b = -5, l = 5)) +     # reset margins with _10_ in every spot
  geom_segment(aes(x = 3.5, y = 0, xend = 3.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 6.5, y = 0, xend = 6.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 9.5, y = 0, xend = 9.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 12.5, y = 0, xend = 12.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 15.5, y = 0, xend = 15.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 18.5, y = 0, xend = 18.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 21.5, y = 0, xend = 21.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 24.5, y = 0, xend = 24.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 27.5, y = 0, xend = 27.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 30.5, y = 0, xend = 30.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 33.5, y = 0, xend = 33.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 36.5, y = 0, xend = 36.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 39.5, y = 0, xend = 39.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  geom_segment(aes(x = 42.5, y = 0, xend = 42.5, yend = 100), linewidth = 0.4, color = "black") +   # demarcate Tx grps
  facet_wrap(~host.ptg, nrow = 2)  
# we see the within-sample DVG community trends. 
#     Firstly, there appears to very little flex in the DVG sub-community relative abundances. Although the Tx's clearly affects the counts (numbered in black), the proportions across treatments remains stable
#     Secondly, within this hyper-stable DVG sub-community, ca09 and tx12 have unique compositions:
#           i. for ca09, the HA and NP DVG segments occur with the highest relative abundance, with PA following as a close third.
#          ii. for tx12, the M and PA DVG segments occur with the highest relative abundance
#          z1. this same trend was observed in our Alpelisib study
#          z2. Leptolyngbya completely 
#     Given the count of DVGs of a particular segment at 18 h.p.i., it may be possible to derive the DVG count of other genomic segments


#_______export the plots...conveniently export plots w/ggsave()
# Set desired save directory
save_directory <- "~/R_gcay_screen"
# save to above directory
ggsave(file.path(save_directory, "TVG.Count.per.Segment.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.1.gcay.data.spread.sgmt, 
       width = 7 , height = 9)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "Mean.TVG.Count.per.Segment.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.2.gcay.data.spread.sgmt, 
       width = 7 , height = 9)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "Mean.DVG.Proportion.per.Segment.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.4.gcay.data.spread.sgmt, 
       width = 7 , height = 9)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "within.sample.prop.dvg.as.pct.of.total.dvg.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       overview.5.gcay.data.spread.sgmt, 
       width = 7 , height = 9)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements





#______________RESULTS 1.1____________per segment level
# KNOWN #1: intimate metabolic cross-talk between flu virus and host metabolism
# KNOWN #2: At 18h.p.i., Alpelisib has significantly different effects on total genomes and DVG proportion at the per segment level
# KNOWN #3: THis effect of Alpelisb (KNOWN#2) also extends to the viral particle level, where it can also have significantly different impacts on TCU titer and the proportion of NCUs
#
# UNKNOWN: the impact of a screen of treatments dosing on flu progeny evolution i.e. SVG vs DVG outcomes AT THE SEGMENT LEVEL
#                   i.   percent change in tvg count for each segment...'count.tvg.sgmt' 
#                  ii.  percent change in relative abundance of dvgs for each segment... 'prop.dvg.sgmt'     
#          i.e. can a treatment induce significantly different effects on the parameters that define DI; 'count.tvg.sgmt' and 'prop.dvg.sgmt'
#               Ha: %change in prop.dvg.sgmt ≠ %change in count.tvg.sgmt
#               H0: %change in prop.dvg.sgmt = %change in count.tvg.sgmt
#
#           IF there is a significant difference, 
#             THEN are the observed changes consistent with the generally accepted definition of DI?
#                 i.  %change in count.tvg.sgmt must be negative in magnitude
#                 ii. %change in prop.dvg may be positive or slightly negative in magnitude




# Order of events
# i. create percent change object...containing 'pct.change.count.tvg' and 'pct.change.prop.dvg'
# ii. plot chart
# iii. hypothesis-test




#________percent change dataframe
# Okay, here's how you do this
# I want to calculate %change in count.tvg between 
#     vehicle controls(DMSO/dH2O) and the other Tx's 
#      but within strains i.e. for ca09 then tx12
#         AND FOR EACH SEGMENT
#	'%change' formula = ((new.value - baseline.value) / baseline.value) * 100
#		where
#		- baseline.value   is count.tvg for DMSO/dH2O
#		- new.value	is count.tvg for 'Tx' ("Alpe" , "MK2206" , "4-OI" , "UK5099" , "Insu" , "Ado" , "Uri" , "Lepto" , "Nosto" , "Oscil" , "Synec" , "Tolyp")
#     i.e. 'count.tvg.sgmt'
#
#
# To do	new.value - baseline.value ...more accurately count.tvg.sgmt - count.tvg.sgmt.baseline
#        I need to create/add a		count.tvg.sgmt.baseline	variable to the df
#	       but populate it w/the the DMSO count.tvg.sgmt for the DMSO-vehicle grps, and w/dH2O count.tvg.sgmt for the dH2O-vehicle grps i.e. DMSO for PB2 PB1 PA HA NP NA M NS ,  and dH2O for PB2 PB1 PA HA NP NA M NS 
#	       the challenge is that ea bioreplicate batch has its own unique count.tvg.baseline of DMSO and dH2O
#           so I need to create a df for each bioreplicate(3x) for each strain(2x) for each vehicle(2x) and for each segment(8x)
#           so     48x dataframes ea for ca09 and tx12...
#		I will do this w/ a combo of    	mutate(...ifelse)     and      fill()            shown in a bit
#
# then bind all 96x df's w/       bind_rows(df1, df2, df3, ... , df96)      48ea for ca09 + tx12
# then plot with ggplot()
# then hypothesis test


# now, make the 48x df's...they'll be bound after into the a 'pct.change' df
#	Instead of 48x df's, I can use the right combo of filtering and grouping to do it in a single df
# REMINDER:
#     these df's are to calculate percent change relative to the vehicle ctrl


# pct.change.genomes.alpvflu1.dmso
# pct.change.genomes.alpvflu1.dh2o


### As I create the datafame,
#       I will assign the unique bioreplicate vehicle-treated (0uM) infection like so:
#
#
# the	   mutate(...ifelse) and fill()	     combo of codes below is to fill-in the new count.tvg.baseline variable for all observations
#	mutate() creates the new variable
#		ifelse() says "if the value in the Alpe.uM.fctr variable is '0' i.e. vehicle-treated, then populate with that observation's count.tvg value
#		if not, populate w/ 'NA'
#	then the fill() cmd fills-in the 'NA' values with the last computed value, which the count.tvg of the '0uM' test group

library(tidyverse)

#____the dmso vehicle object___
pct.change.genomes.alpvflu1.dmso <-
  dvg.svg.bind %>%
  filter(host.ptg %in% c("MDCK.CA09(H1N1)", "MDCK.TX12(H3N2)" ),
         veh.Tx=="DMSO"
  ) %>%
  group_by(host.ptg, bioreplicate, Tx, geno.sgmt) %>%
  summarise(count.tvg.sgmt = n(),
            count.dvg.sgmt = sum(geno.class == "DVG"),
            prop.dvg.sgmt = (sum(geno.class == "DVG")) / n()    ) %>%
  ungroup() %>%
  
  mutate(count.tvg.sgmt.baseline = ifelse(Tx =="DMSO" , count.tvg.sgmt , NA ) ) %>%           # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.sgmt.baseline, .direction = "down") %>%                                                 # fills-in the 'NA' values with the last computed value, which the count.tvg.sgmt of the 'DMSO' test group
  mutate(pct.change.count.tvg.sgmt = ((count.tvg.sgmt - count.tvg.sgmt.baseline)/count.tvg.sgmt.baseline) * 100) %>%         # create a variable to hold the percent change of TVG count between the vehicle-treated group and each Tx
  
  mutate(prop.dvg.sgmt.baseline = ifelse(Tx =="DMSO" , prop.dvg.sgmt , NA ) ) %>%           # create column for baseline DVG proportion; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(prop.dvg.sgmt.baseline, .direction = "down") %>%                                                 # fills-in the 'NA' values with the last computed value, which the prop.dvg.sgmt of the 'DMSO' test group
  mutate(pct.change.prop.dvg.sgmt = ((prop.dvg.sgmt - prop.dvg.sgmt.baseline)/prop.dvg.sgmt.baseline) * 100)        # create a variable to hold the percent change of DVG proportion between the vehicle-treated group and each Tx

# verify the number of records in the 'pct.change' object
# 2x strains * 3x bioreplicates * 8x segments * 6x treatments w/DMSO vehicle = 288
#   there should be 288 records in    View(pct.change.genomes.alpvflu1.dmso )

#____the dh2o vehicle object___
pct.change.genomes.alpvflu1.dh2o <-
  dvg.svg.bind %>%
  filter(host.ptg %in% c("MDCK.CA09(H1N1)", "MDCK.TX12(H3N2)" ),
         veh.Tx=="dH2O"
  ) %>%
  group_by(host.ptg, bioreplicate, Tx, geno.sgmt) %>%
  summarise(count.tvg.sgmt = n(),
            count.dvg.sgmt = sum(geno.class == "DVG"),
            prop.dvg.sgmt = (sum(geno.class == "DVG")) / n()    ) %>%
  ungroup() %>%
  
  mutate(count.tvg.sgmt.baseline = ifelse(Tx =="dH2O" , count.tvg.sgmt , NA ) ) %>%           # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.sgmt.baseline, .direction = "down") %>%                                                 # fills-in the 'NA' values with the last computed value, which the count.tvg.sgmt of the 'dH2O' test group
  mutate(pct.change.count.tvg.sgmt = ((count.tvg.sgmt - count.tvg.sgmt.baseline)/count.tvg.sgmt.baseline) * 100) %>%         # create a variable to hold the percent change of TVG count between the vehicle-treated group and each Tx
  
  mutate(prop.dvg.sgmt.baseline = ifelse(Tx =="dH2O" , prop.dvg.sgmt , NA ) ) %>%           # create column for baseline DVG proportion; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(prop.dvg.sgmt.baseline, .direction = "down") %>%                                                 # fills-in the 'NA' values with the last computed value, which the prop.dvg.sgmt of the 'dH2O' test group
  mutate(pct.change.prop.dvg.sgmt = ((prop.dvg.sgmt - prop.dvg.sgmt.baseline)/prop.dvg.sgmt.baseline) * 100)         # create a variable to hold the percent change of DVG proportion between the vehicle-treated group and each Tx

# verify the number of records in the 'pct.change' object
# 2x strains * 3x bioreplicates * 8x segments * 9x treatments w/dH2O vehicle = 432
#   there should be 432 records in    View(pct.change.genomes.alpvflu1.dh2o )
#     but there are 429 ... because some segments came back with 0 
#        TX12 "Uri" bioreplicate #1 had had zero recovered HA and NA segments
#        TX12 "Uri" bioreplicate #2 had had zero recovered NA segments
#        429 + 3 = 432   so  ALL CLEAR


#_____bind %change data_____
pct.change.genomes.sgmt.alpvflu1.bind <-
  bind_rows(pct.change.genomes.alpvflu1.dmso,
            pct.change.genomes.alpvflu1.dh2o) 

#____________________barchart________________________________
# First we plot point charts--NOT barcharts--
#    This is to get a sense of heterogeneity following normalization via %change calculation
#         within bioreplicates

# pct.change.count.tvg.sgmt and pct.change.prop.dvg.sgmt
#   the per segment charts are pretty packed, so I'm splitting up the segments

#___ca09 plot
chart.pct.change.sgmt.tvg_v_prop.dvg_ca09 <-
  pct.change.genomes.sgmt.alpvflu1.bind %>%
  filter(host.ptg=="MDCK.CA09(H1N1)") %>%
  ggplot( aes(x = Tx, shape = bioreplicate)) +
  geom_point(aes(y = pct.change.count.tvg.sgmt, color = "pct.change.count.tvg.sgmt", shape= bioreplicate), position = position_jitter(width = 0.05), size = 3, alpha = 0.5) +
  geom_point(aes(y = pct.change.prop.dvg.sgmt, color = "pct.change.prop.dvg.sgmt", shape= bioreplicate), position = position_jitter(width = 0.05), size = 3, alpha = 0.5) +
  labs(y = "% Change", x= "Treatment", title = "MDCK.CA09(H1N1)") +
  geom_text(aes(y = pct.change.count.tvg.sgmt, label = round(pct.change.count.tvg.sgmt, 1)), vjust = -0.5, size = 2.5, nudge_x = -0.2) +     # an alternate to nudge_x is      position = position_jitter(width = 0.4)
  geom_text(aes(y = pct.change.prop.dvg.sgmt, label = round(pct.change.prop.dvg.sgmt, 1)), vjust = -0.5, size = 2.5, nudge_x = 0.2) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7)) +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "DI Marker", direction = "horizontal"),         # set 'title =NULL' to delete legend title...set 'direction=horizontal' to flatten out legend
         shape = guide_legend(title = "Bioreplicate", direction = "horizontal")    ) +
  scale_color_manual(values = c("pct.change.count.tvg.sgmt" = "blue", "pct.change.prop.dvg.sgmt" = "red")) +
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 18)) +
  facet_wrap(~geno.sgmt,
             scales = "free_y"         # FYI: you can't "free_y" and still set ylim() ... silence this line to test out the ylim() codes below
  ) 
  #ylim(-100, 8000)          # get a better look at PB2 PB1 PA HA and NP segments
  #ylim(-100, 300)          # get a better look at NA, M, and NS segments

#___tx12 plot
chart.pct.change.sgmt.tvg_v_prop.dvg_tx12 <-
  pct.change.genomes.sgmt.alpvflu1.bind %>%
  filter(host.ptg=="MDCK.TX12(H3N2)") %>%
  ggplot( aes(x = Tx, shape = bioreplicate)) +
  geom_point(aes(y = pct.change.count.tvg.sgmt, color = "pct.change.count.tvg.sgmt", shape= bioreplicate), position = position_jitter(width = 0.05), size = 3, alpha = 0.5) +
  geom_point(aes(y = pct.change.prop.dvg.sgmt, color = "pct.change.prop.dvg.sgmt", shape= bioreplicate), position = position_jitter(width = 0.05), size = 3, alpha = 0.5) +
  labs(y = "% Change", x= "Treatment", title = "MDCK.TX12(H3N2)") +
  geom_text(aes(y = pct.change.count.tvg.sgmt, label = round(pct.change.count.tvg.sgmt, 1)), vjust = -0.5, size = 2.5, nudge_x = -0.2) +     # an alternate to nudge_x is      position = position_jitter(width = 0.4)
  geom_text(aes(y = pct.change.prop.dvg.sgmt, label = round(pct.change.prop.dvg.sgmt, 1)), vjust = -0.5, size = 2.5, nudge_x = 0.2) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7)) +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10)) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "DI Marker", direction = "horizontal"),         # set 'title =NULL' to delete legend title...set 'direction=horizontal' to flatten out legend
         shape = guide_legend(title = "Bioreplicate", direction = "horizontal")    ) +
  scale_color_manual(values = c("pct.change.count.tvg.sgmt" = "blue", "pct.change.prop.dvg.sgmt" = "red")) +
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 18)) +
  facet_wrap(~geno.sgmt,
             #scales = "free_y"         # FYI: you can't "free_y" and still set ylim() ... silence this line to test out the ylim() codes below
  ) 
  #ylim(-100, 1200)          # get a better look at PB2 PB1 PA and M segments
  #ylim(-100, 400)          # get a better look at the other segments




#_______export the plot...conveniently export plots w/ggsave()
# Set desired save directory
save_directory <- "~/R_gcay_screen"
# save to above directory
ggsave(file.path(save_directory, "Percent.Change.per.Segment.tvg_v_prop.dvg_ca09.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       chart.pct.change.sgmt.tvg_v_prop.dvg_ca09, 
       width = 7 , height = 6)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "Percent.Change.per.Segment.tvg_v_prop.dvg_tx12.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       chart.pct.change.sgmt.tvg_v_prop.dvg_tx12, 
       width = 7 , height = 6)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements


#____________bioreplicate-averaged barchart____________
# mean pct.change.count.tvg.sgmt and mean pct.change.prop.ncu.sgmt
#
#
# create the mean-computed object first, then call it with ggplot()
pct.change.genomes.tvg_v_prop.dvg.mean.sgmt <-
  pct.change.genomes.sgmt.alpvflu1.bind %>% 
  group_by(host.ptg, geno.sgmt, Tx   # grouping for mean of bioreplicates
  ) %>% 
  summarise(n = n(),
            mean.pct.tvg.sgmt = mean(pct.change.count.tvg.sgmt),     # mean TVG count (dosing experiment bioreplicates)
            sd.tvg.sgmt = sd(pct.change.count.tvg.sgmt),               # sd
            sem.tvg.sgmt = sd.tvg.sgmt /sqrt(n()),          #sem
            
            mean.prop.dvg.sgmt = mean(pct.change.prop.dvg.sgmt),     # mean DVG relative abundance (dosing experiment bioreplicates)
            sd.prop.dvg.sgmt = sd(pct.change.prop.dvg.sgmt),           #sd
            sem.prop.dvg.sgmt = sd.prop.dvg.sgmt /sqrt(n())   #sem
  ) %>%
  ungroup() %>%
  rename(
    pct.change.count.tvg.sgmt = mean.pct.tvg.sgmt,                                           # rename variable 
    pct.change.prop.dvg.sgmt = mean.prop.dvg.sgmt,                                           # rename variable
  )


#______plot barchart
#chart.pct.change.tvg_v_prop.dvg.mean.sgmt <-
pct.change.genomes.tvg_v_prop.dvg.mean.sgmt %>%
  filter(!Tx %in% c("DMSO", "dH2O")) %>%         # the plot is really busy, so I'm dropping the vehicle controls
  ggplot(aes(x = Tx)) +
  geom_bar(aes(y = pct.change.count.tvg.sgmt, fill = "pct.change.count.tvg.sgmt"), 
           position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), 
           stat = "identity", alpha = 0.4) +
  geom_errorbar(aes(y = pct.change.count.tvg.sgmt, 
                    ymin = pct.change.count.tvg.sgmt - sem.tvg.sgmt, 
                    ymax = pct.change.count.tvg.sgmt + sem.tvg.sgmt),
                position = position_dodge(width = 0.9), width = 0.2, alpha = 0.5) +
  geom_bar(aes(y = pct.change.prop.dvg.sgmt, fill = "pct.change.prop.dvg.sgmt"), 
           position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), 
           stat = "identity", alpha = 0.5) +
  geom_errorbar(aes(y = pct.change.prop.dvg.sgmt, 
                    ymin = pct.change.prop.dvg.sgmt - sem.prop.dvg.sgmt, 
                    ymax = pct.change.prop.dvg.sgmt + sem.prop.dvg.sgmt),
                position = position_dodge(width = 0.9), width = 0.2, alpha = 0.5) +
  labs(y = "% Change", x= "Alpelisib (uM)", fill = NULL) +  # Remove legend title with 'fill=NULL'...it's clunky and the 
  geom_text(aes(y = pct.change.count.tvg.sgmt, label = round(pct.change.count.tvg.sgmt, 0)), vjust = -0.5, size = 2, nudge_x = -0.2) +     # an alternate to nudge_x is      position = position_jitter(width = 0.4)
  geom_text(aes(y = pct.change.prop.dvg.sgmt, label = round(pct.change.prop.dvg.sgmt, 0)), vjust = -0.5, size = 2, nudge_x = 0.2) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7, angle=45)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("pct.change.count.tvg.sgmt" = "blue", "pct.change.prop.dvg.sgmt" = "red")) +
  facet_grid(geno.sgmt~host.ptg,
             scales = "free_y")         # the base plot shows that %change of prop.dvg gets really high, the plot need to be fragged out
#coord_cartesian(ylim = c(-100, 100))    # there's no one-size-fit-all for this plot ... "free_y" is the best bet ... if you start to frag out the plots, it becomes a mess ...just rely on the asterisks to plainly communicate where the differences were significant


#
#
# don't yet export the plot...you may need to add asterisks for significance
# NEXT do the hypothesis testing to determine significance,
#    then add the asterisks and line segments denoting significance to the chart object
#
# promising treatments
#   Normally this section would go something like this:
#       DON'T waste time testing every pair...that comes out to 100+ tests for the per-segment data
#       focus on promising test groups
#       these are those treatments where in a point chart, there is NO overlap between values of both groups
#         don't rely on SD or SEM error bars to determine overlap. 
#             SD and SEM can show no overlap even when some of the individual data points overlap
# 
#   In the drug screen however, a good bit look promising so just do them all
#      automate as much as you can
#
#



#___________________
#____________Hypothesis Testing________
#
#
#
# Define the values for geno.sgmt
geno_values <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

# Initialize an empty list to store the results
#   create a list for each Tx for each strain
#   REMINDER: testing is on percent.change values
#    so vehicle treatments (DMSO, dH2O) which serve as the baseline value have 0% change
#      0% change n=3 times means they have no distribution and are not testable
#       so testing is just within the non vehicle Tx groups
#        Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
if(FALSE) {
  
  result_list <- list()   
  
}                            # if(FALSE) ends here

result_list.dmso.ca09 <- list()
result_list.dh2o.ca09 <- list()
result_list.insu.ca09 <- list()
result_list.alpe.ca09 <- list()
result_list.mk2206.ca09 <- list()
result_list.4oi.ca09 <- list()
result_list.uk5099.ca09 <- list()
result_list.ado.ca09 <- list()
result_list.uri.ca09 <- list()
result_list.favp.ca09 <- list()
result_list.lepto.ca09 <- list()
result_list.nosto.ca09 <- list()
result_list.oscil.ca09 <- list()
result_list.synec.ca09 <- list()
result_list.tolyp.ca09 <- list()

result_list.dmso.tx12 <- list()
result_list.dh2o.tx12 <- list()
result_list.insu.tx12 <- list()
result_list.alpe.tx12 <- list()
result_list.mk2206.tx12 <- list()
result_list.4oi.tx12 <- list()
result_list.uk5099.tx12 <- list()
result_list.ado.tx12 <- list()
result_list.uri.tx12 <- list()
result_list.favp.tx12 <- list()
result_list.lepto.tx12 <- list()
result_list.nosto.tx12 <- list()
result_list.oscil.tx12 <- list()
result_list.synec.tx12 <- list()
result_list.tolyp.tx12 <- list()


# Use below loop to manually populate the empty lists generated above
# Loop through each value of geno.sgmt
#   !  ALERT    ALERT     ALERT    !
#        I had to highlight the entire _for_ loop before it ran correctly...just using Ctrl-Enter didn't work
#        ALSO, be sure to use the percent change object, not the base values object

# sample code ... DO NOT RUN
#   if(FALSE) {} is so a script run (Ctrl-Shift-S) will bypass it this sample code
if(FALSE) {
  
  for (geno_value in geno_values) {
    # Construct the filter condition dynamically
    filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
    
    # Run the command and store the result in the empty list created earlier
    result_list.40.tx12[[geno_value]] <- gcay.alpe.count.sgmt %>%
      filter(host.ptg == "MDCK.TX12(H3N2)", Alpe.uM.fctr == "40", !!rlang::parse_expr(filter_condition)) %>%
      select(count.tvg.sgmt, prop.dvg.sgmt)
  }
  
}     # if(FALSE) ends here

# The output is a list which contains tibbles (data frames) for each segments. 
# Each tibble contains two columns: pct.change.count.tvg.sgmt and pct.change.prop.dvg.sgmt.
# To access data within a specific segment e.g. "PB2" ,  use result_list.uri.tx12$PB2
#    and to access data from a two-column tibble list, use double square brackets. For example, to access the results for the first column of the "PB2" segment, you can use result_list.uri.tx12$PB2[[1]].

# result_list.dmso.ca09
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.dmso.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "DMSO", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


# result_list.dh2o.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.dh2o.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "dH2O", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.insu.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.insu.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Insu", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.alpe.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.alpe.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Alpe", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.mk2206.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.mk2206.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "MK2206", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.4oi.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.4oi.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "4-OI", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.uk5099.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.uk5099.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "UK5099", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.ado.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.ado.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Ado", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.uri.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.uri.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Uri", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.favp.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.favp.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Favp", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.lepto.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.lepto.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Lepto", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.nosto.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.nosto.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Nosto", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.oscil.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.oscil.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Oscil", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.synec.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.synec.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Synec", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.tolyp.ca09 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.tolyp.ca09[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.CA09(H1N1)", Tx == "Tolyp", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}

#_______________________________
#result_list.dmso.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.dmso.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "DMSO", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.dh2o.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.dh2o.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "dH2O", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.insu.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.insu.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Insu", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.alpe.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.alpe.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Alpe", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.mk2206.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.mk2206.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "MK2206", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.4oi.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.4oi.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "4-OI", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.uk5099.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.uk5099.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "UK5099", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.ado.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.ado.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Ado", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.uri.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.uri.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Uri", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.favp.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.favp.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Favp", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.lepto.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.lepto.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Lepto", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.nosto.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.nosto.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Nosto", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.oscil.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.oscil.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Oscil", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.synec.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.synec.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Synec", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}


#result_list.tolyp.tx12 <- list()
for (geno_value in geno_values) {
  # Construct the filter condition dynamically
  filter_condition <- paste("geno.sgmt == '", geno_value, "'", sep = "")
  
  # Run the command and store the result in the empty list created earlier
  result_list.tolyp.tx12[[geno_value]] <- pct.change.genomes.sgmt.alpvflu1.bind %>%
    filter(host.ptg == "MDCK.TX12(H3N2)", Tx == "Tolyp", !!rlang::parse_expr(filter_condition)) %>%
    select(pct.change.count.tvg.sgmt, pct.change.prop.dvg.sgmt)
}



#___normality
# Apply shapiro.test() to each tibble in result_list
#   The variables are percent.change, so vehicle groups (dH2O  DMSO) will be zero ... don't test them
#    Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
#
# There are different codes to run normality test for the 1- or 2-column input objects
#     Option 1: run normality test on a list of 1-column tibbles, and output the result
#     Option 2: run normality test on a list of 2-column tibbles, and output the results
# Sample codes for both options are below
#    don't run these sample codes ... just copy and paste as needed
#    if(FALSE) {} is wrapped around the sample codes so that a running script (Ctrl-Shift-S) will bypass these sample codes
# ! ALERT    ALERT     ALERT   !
#        because of the "tibbles in a list" format of the input files, 
#         the entire operation will fail if even a single tibble has an issue
#           the two main issues are less than n=3 replicates, or bioreplicates having the same value
#            the codes were written output 'NA' for those segments in the list that violate any of 
#             these requirements so that the full operation can still run and output results for the other compliant segments
#
#____shapiro.test CODES EXPLAINED______________
# We use lapply() to iterate over each tibble in result_list.uri.tx12.
# For each tibble, we use sapply() to iterate over each column.
# Within the inner sapply() loop, we perform the Shapiro-Wilk test on each column.
# 
# We include logic to handle cases where all values are identical:
#   We use length(unique(col)) > 1 to check if there is more than one unique value in the column. If all values are identical, length(unique(col)) will be 1.
# If all values are identical, we issue a warning and assign NA to the corresponding result in the output list.
# Otherwise, we perform the Shapiro-Wilk test on the column
# 
# We include logic to handle cases where the sample size is not within the required range:
#   We're using nrow(tbl) to determine the number of rows (sample size) in each tibble.
#     	We're checking if the sample size is within the required range (3 to 5000) before performing the Shapiro-Wilk test.
# If the sample size is not within the required range, we issue a warning and assign NA to the corresponding result in the output list.
# 
# The output is a nested list where each element corresponds to a tibble, and each nested element corresponds to a column in that tibble. The nested elements contain the results of the Shapiro-Wilk test for each column.
# 
# Access the Shapiro-Wilk test results as follows:
#    from the single column tibble list ... use '$' ... e.g. for "PB2"   use result_list.uri.tx12$PB2.
#    from the two-column tibble list ... use '$' twice ... for "PB2" prop.dvg.sgmt column, use result_list.uri.tx12$PB2$pct.change.prop.dvg.sgmt

#    if(FALSE) {} is wrapped around the sample codes so that a running script (Ctrl-Shift-S) will bypass these sample codes

if(FALSE) {    #if(FALSE) starts here
  
  # run normality test and output the results of both columns in the list of tibbles
  # "equal values" ERROR fixed
  # "minimum replicates" ERROR fixed
  #
  # tx12 lepto test run
  shapiro_tests.lepto.tx12 <- lapply(result_list.lepto.tx12, function(tbl) {
    sapply(tbl, function(col) {
      if (length(unique(col)) > 1) { # Check if there is more than one unique value
        if (length(col) >= 3 && length(col) <= 5000) {
          shapiro.test(col)
        } else {
          warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
          NA
        }
      } else {
        warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
        NA
      }
    })
  })
  
  
  # The below code is similar to above code but for a single column at a time
  #   normality test and output the result of a single column in the list of tibbles
  #
  # Define the column name you want to test
  column_name <- "pct.change.count.tvg.sgmt"
  
  # Apply Shapiro-Wilk test to the specified column of each tibble in result_list.lepto.tx12
  shapiro_tests.column-name.lepto.tx12 <- lapply(result_list.lepto.tx12, function(tbl) {
    if (column_name %in% names(tbl)) {
      col <- tbl[[column_name]]
      if (length(unique(col)) > 1) { # Check if there is more than one unique value
        if (length(col) >= 3 && length(col) <= 5000) {
          shapiro.test(col)
        } else {
          warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", column_name)
          NA
        }
      } else {
        warning("All 'x' values are identical for column ", column_name, ". Shapiro-Wilk test cannot be performed.")
        NA
      }
    } else {
      warning("Column ", column_name, " does not exist in the tibble.")
      NA
    }
  })
  
}      # if(FALSE) ends here


# We proceed on the 2-column tibble route
# REMINDER: test groups are 
#  Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
#ca09

shapiro_tests.insu.ca09 <- lapply(result_list.insu.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.alpe.ca09 <- lapply(result_list.alpe.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.mk2206.ca09 <- lapply(result_list.mk2206.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.4oi.ca09 <- lapply(result_list.4oi.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.uk5099.ca09 <- lapply(result_list.uk5099.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.ado.ca09 <- lapply(result_list.ado.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.uri.ca09 <- lapply(result_list.uri.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.favp.ca09 <- lapply(result_list.favp.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.lepto.ca09 <- lapply(result_list.lepto.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.nosto.ca09 <- lapply(result_list.nosto.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.oscil.ca09 <- lapply(result_list.oscil.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.synec.ca09 <- lapply(result_list.synec.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.tolyp.ca09 <- lapply(result_list.tolyp.ca09, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})

#
# RESULTS:  which groups are normality-negative ... they will undergo non-parametric testing
#           all data distribution results are after tx12 normality testing
#
#tx12

shapiro_tests.insu.tx12 <- lapply(result_list.insu.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.alpe.tx12 <- lapply(result_list.alpe.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.mk2206.tx12 <- lapply(result_list.mk2206.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.4oi.tx12 <- lapply(result_list.4oi.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.uk5099.tx12 <- lapply(result_list.uk5099.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.ado.tx12 <- lapply(result_list.ado.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.uri.tx12 <- lapply(result_list.uri.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.favp.tx12 <- lapply(result_list.favp.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.lepto.tx12 <- lapply(result_list.lepto.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.nosto.tx12 <- lapply(result_list.nosto.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.oscil.tx12 <- lapply(result_list.oscil.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.synec.tx12 <- lapply(result_list.synec.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})


shapiro_tests.tolyp.tx12 <- lapply(result_list.tolyp.tx12, function(tbl) {
  sapply(tbl, function(col) {
    if (length(unique(col)) > 1) { # Check if there is more than one unique value
      if (length(col) >= 3 && length(col) <= 5000) {
        shapiro.test(col)
      } else {
        warning("Sample size is not within the required range (3 to 5000). Shapiro-Wilk test cannot be performed for column ", deparse(substitute(col)))
        NA
      }
    } else {
      warning("All 'x' values are identical for column ", deparse(substitute(col)), ". Shapiro-Wilk test cannot be performed.")
      NA
    }
  })
})

#
# 
# RESULTS:  
# below groups are normality-negative ... they'll undergo non-parametric Wilcoxon Signed-Rank Testing
# the unlisted normality-positive groups will undergo parametric T-Test
# 
#
#
# shapiro_tests.alpe.ca09		pb2 pa ha np
# shapiro_tests.mk2206.ca09	pb1 pa
# shapiro_tests.uk5099.ca09	pb2 pb1
# shapiro_tests.ado.ca09		pb1
# shapiro_tests.uri.ca09		pb2 pa
# shapiro_tests.favp.ca09		ha np na m ns
# shapiro_tests.lepto.ca09	ha
# shapiro_tests.nosto.ca09	na m
# shapiro_tests.oscil.ca09	pb1 np
# shapiro_tests.tolyp.ca09	pa

# shapiro_tests.insu.tx12		pb1 ha
# shapiro_tests.mk2206.tx12	pb1 ha np
# shapiro_tests.4oi.tx12		pb1
# shapiro_tests.uk5099.tx12	pb1
# shapiro_tests.ado.tx12		pa na m
# shapiro_tests.uri.tx12		pa ns
# shapiro_tests.favp.tx12		ha na
# shapiro_tests.lepto.tx12	ns.... NOTE: ha + na are not testable(n=2 DVGs)
# shapiro_tests.nosto.tx12	ha np na m ns
# shapiro_tests.oscil.tx12	pb2 na
# shapiro_tests.tolyp.tx12	pa ha na
#


#__________equality of variances________________
# test equality of variance
# Not necessary here:
#   equal variances are not a required assumption for paired T-Test, only that data is normality-positive
#       normality-negative data will undergo paired sample Wilcoxon Signed Rank Test
#             wilcox.test(Grp1, Grp2,  
#                         paired = TRUE,         # One sample, paired observations
#                         alternative = "two.sided")
#       normality-positive data will undergo paired sample T-Test
#             t.test(Grp1, Grp2,
#                    paired = TRUE,               # One sample, paired observations
#                    var.equal = TRUE,           # unrequired assumption
#                    alternative = "two.sided")
# p=0.3


#_____parametric or non-parametric hypothesis test?__________
#
# normality-negative data sets will undergo Wilcoxon Signed Rank Test
# normality-positive, variance-unequal data sets will undergo Welch's T-Test...IF unpaired...not applicable for Paired T-Test
# normality-positive, variance-equal data sets will undergo Student's T-Test
# ALL tests are paired, since test variables are from the same sample

#_______Paired one-sample two-variable ttest
# REMINDER: test groups are 
#  Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp


#___ca_insu___
# normality negative
#   NONE

# normality positive
#  pb2 pb1 pa ha np na m ns
#ttest.pct.change.genomes.pb2.insu.ca09 <- 
t.test(result_list.insu.ca09$PB2$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$PB2$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.25

#ttest.pct.change.genomes.pb1.insu.ca09 <- 
t.test(result_list.insu.ca09$PB1$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$PB1$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.21

#ttest.pct.change.genomes.pa.insu.ca09 <- 
t.test(result_list.insu.ca09$PA$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$PA$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.22

#ttest.pct.change.genomes.ha.insu.ca09 <- 
t.test(result_list.insu.ca09$HA$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$HA$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.54

#ttest.pct.change.genomes.np.insu.ca09 <- 
t.test(result_list.insu.ca09$NP$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$NP$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.01 *******************************

#ttest.pct.change.genomes.na.insu.ca09 <- 
t.test(result_list.insu.ca09$`NA`$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$`NA`$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.006********************

#ttest.pct.change.genomes.m.insu.ca09 <- 
t.test(result_list.insu.ca09$M$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$M$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.8

#ttest.pct.change.genomes.ns.insu.ca09 <- 
t.test(result_list.insu.ca09$NS$pct.change.prop.dvg.sgmt, 
       result_list.insu.ca09$NS$pct.change.count.tvg.sgmt,
       paired = TRUE,               # One sample, paired observations
       #var.equal = TRUE,           # unrequired assumption
       alternative = "two.sided")
# p=0.29

#___ca_alpe___
# normality-negative
# pb2 pa ha np
#
# normality-positive
# pb1 na m ns

# ttest.pct.change.genomes.pb1.alpe.ca09 <- 
t.test(result_list.alpe.ca09$PB1$pct.change.prop.dvg.sgmt,
            result_list.alpe.ca09$PB1$pct.change.count.tvg.sgmt,  
            paired = TRUE, 
            alternative = "two.sided")
#p= 0.32

# ttest.pct.change.genomes.na.alpe.ca09 <- 
t.test(result_list.alpe.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.alpe.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.08

# ttest.pct.change.genomes.m.alpe.ca09 <- 
t.test(result_list.alpe.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.alpe.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.2

# ttest.pct.change.genomes.ns.alpe.ca09 <- 
t.test(result_list.alpe.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.alpe.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.14

#___ca_mk2206___
# normality-negative
# pb1 pa
#
# normality-positive
# pb2 ha np na m ns

# ttest.pct.change.genomes.pb2.mk2206.ca09 <- 
t.test(result_list.mk2206.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.mk2206.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.2

# ttest.pct.change.genomes.ha.mk2206.ca09 <- 
t.test(result_list.mk2206.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.mk2206.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.07

# ttest.pct.change.genomes.np.mk2206.ca09 <- 
t.test(result_list.mk2206.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.mk2206.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.053

# ttest.pct.change.genomes.na.mk2206.ca09 <- 
t.test(result_list.mk2206.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.mk2206.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.32

# ttest.pct.change.genomes.m.mk2206.ca09 <- 
t.test(result_list.mk2206.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.mk2206.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.02 ******************************

# ttest.pct.change.genomes.ns.mk2206.ca09 <- 
t.test(result_list.mk2206.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.mk2206.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.3

#___ca_4oi___
# normality-negative
# NONE
#
# normality-positive
# pb2 pb1 pa ha np na m ns

# ttest.pct.change.genomes.pb2.4oi.ca09 <- 
t.test(result_list.4oi.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.07

# ttest.pct.change.genomes.pb1.4oi.ca09 <- 
t.test(result_list.4oi.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.10

# ttest.pct.change.genomes.pa.4oi.ca09 <- 
t.test(result_list.4oi.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.042 *******************************

# ttest.pct.change.genomes.ha.4oi.ca09 <- 
t.test(result_list.4oi.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.08

# ttest.pct.change.genomes.np.4oi.ca09 <- 
t.test(result_list.4oi.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.052

# ttest.pct.change.genomes.na.4oi.ca09 <- 
t.test(result_list.4oi.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.11

# ttest.pct.change.genomes.m.4oi.ca09 <- 
t.test(result_list.4oi.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.027 ***************

# ttest.pct.change.genomes.ns.4oi.ca09 <- 
t.test(result_list.4oi.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.4oi.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.014 *******************


#___ca_uk5099___
# normality-negative
# pb2 pb1
#
# normality-positive
# pa ha np na m ns

# ttest.pct.change.genomes.pa.uk5099.ca09 <- 
t.test(result_list.uk5099.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.uk5099.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.14

# ttest.pct.change.genomes.ha.uk5099.ca09 <- 
t.test(result_list.uk5099.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.uk5099.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 7.191e-05  0.000071 ********************

# ttest.pct.change.genomes.np.uk5099.ca09 <- 
t.test(result_list.uk5099.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.uk5099.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.014 *********************

# ttest.pct.change.genomes.na.uk5099.ca09 <- 
t.test(result_list.uk5099.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.uk5099.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.011 **********************

# ttest.pct.change.genomes.m.uk5099.ca09 <- 
t.test(result_list.uk5099.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.uk5099.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.19

# ttest.pct.change.genomes.ns.uk5099.ca09 <- 
t.test(result_list.uk5099.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.uk5099.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.055

#___ca_ado___
# normality-negative
# pb1
#
# normality-positive
# pb2 pa ha np na m ns

# ttest.pct.change.genomes.pb2.ado.ca09 <- 
t.test(result_list.ado.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.14

# ttest.pct.change.genomes.pa.ado.ca09 <- 
t.test(result_list.ado.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.15

# ttest.pct.change.genomes.ha.ado.ca09 <- 
t.test(result_list.ado.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.25

# ttest.pct.change.genomes.np.ado.ca09 <- 
t.test(result_list.ado.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.20

# ttest.pct.change.genomes.na.ado.ca09 <- 
t.test(result_list.ado.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.12

# ttest.pct.change.genomes.m.ado.ca09 <- 
t.test(result_list.ado.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.67

# ttest.pct.change.genomes.ns.ado.ca09 <- 
t.test(result_list.ado.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.ado.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.24

#___ca_uri___
# normality-negative
# pb2 pa
#
# normality-positive
# pb1 ha np na m ns

# ttest.pct.change.genomes.pb1.uri.ca09 <- 
t.test(result_list.uri.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.uri.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.24

# ttest.pct.change.genomes.ha.uri.ca09 <- 
t.test(result_list.uri.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.uri.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.044 ******************************

# ttest.pct.change.genomes.np.uri.ca09 <- 
t.test(result_list.uri.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.uri.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.023 *********************

# ttest.pct.change.genomes.na.uri.ca09 <- 
t.test(result_list.uri.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.uri.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.003 **************************8

# ttest.pct.change.genomes.m.uri.ca09 <- 
t.test(result_list.uri.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.uri.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.99

# ttest.pct.change.genomes.ns.uri.ca09 <- 
t.test(result_list.uri.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.uri.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.13

#___ca_favp___
# normality-negative
# ha np na m ns
#
# normality-positive
# pb2 pb1 pa

# ttest.pct.change.genomes.pb2.favp.ca09 <- 
t.test(result_list.favp.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.favp.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0430   ***************************

# ttest.pct.change.genomes.pb1.favp.ca09 <- 
t.test(result_list.favp.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.favp.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0438  ***************************

# ttest.pct.change.genomes.pa.favp.ca09 <- 
t.test(result_list.favp.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.favp.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.025 ***************************

#___ca_lepto___
# normality-negative
# ha
#
# normality-positive
# pb2 pb1 pa np na m ns

# ttest.pct.change.genomes.pb2.lepto.ca09 <- 
t.test(result_list.lepto.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.036 ************************

# ttest.pct.change.genomes.pb1.lepto.ca09 <- 
t.test(result_list.lepto.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.17

# ttest.pct.change.genomes.pa.lepto.ca09 <- 
t.test(result_list.lepto.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0017 ***********************

# ttest.pct.change.genomes.np.lepto.ca09 <- 
t.test(result_list.lepto.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.040 *********************

# ttest.pct.change.genomes.na.lepto.ca09 <- 
t.test(result_list.lepto.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0054  *********************

# ttest.pct.change.genomes.m.lepto.ca09 <- 
t.test(result_list.lepto.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0082  *********************

# ttest.pct.change.genomes.ns.lepto.ca09 <- 
t.test(result_list.lepto.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.lepto.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.012  *********************

#___ca_nosto___
# normality-negative
# na m
#
# normality-positive
# pb2 pb1 pa ha np ns

# ttest.pct.change.genomes.pb2.nosto.ca09 <- 
t.test(result_list.nosto.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.nosto.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.052

# ttest.pct.change.genomes.pb1.nosto.ca09 <- 
t.test(result_list.nosto.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.nosto.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0015 *********************

# ttest.pct.change.genomes.pa.nosto.ca09 <- 
t.test(result_list.nosto.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.nosto.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.025 ****************************

# ttest.pct.change.genomes.ha.nosto.ca09 <- 
t.test(result_list.nosto.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.nosto.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.024 ************************

# ttest.pct.change.genomes.np.nosto.ca09 <- 
t.test(result_list.nosto.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.nosto.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.011 ****************************

# ttest.pct.change.genomes.ns.nosto.ca09 <- 
t.test(result_list.nosto.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.nosto.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.071

#___ca_oscil___
# normality-negative
# pb1 np
#
# normality-positive
# pb2 pa ha na m ns

# ttest.pct.change.genomes.pb2.oscil.ca09 <- 
t.test(result_list.oscil.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.oscil.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.087

# ttest.pct.change.genomes.pa.oscil.ca09 <- 
t.test(result_list.oscil.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.oscil.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.003 *****************************

# ttest.pct.change.genomes.ha.oscil.ca09 <- 
t.test(result_list.oscil.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.oscil.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.001 **************************

# ttest.pct.change.genomes.na.oscil.ca09 <- 
t.test(result_list.oscil.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.oscil.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.078

# ttest.pct.change.genomes.m.oscil.ca09 <- 
t.test(result_list.oscil.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.oscil.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.57

# ttest.pct.change.genomes.ns.oscil.ca09 <- 
t.test(result_list.oscil.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.oscil.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.98

#___ca_synec___
# normality-negative
# NONE
#
# normality-positive
# pb2 pb1 pa ha np na m ns

# ttest.pct.change.genomes.pb2.synec.ca09 <- 
t.test(result_list.synec.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.11

# ttest.pct.change.genomes.pb1.synec.ca09 <- 
t.test(result_list.synec.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.12

# ttest.pct.change.genomes.pa.synec.ca09 <- 
t.test(result_list.synec.ca09$PA$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.079

# ttest.pct.change.genomes.ha.synec.ca09 <- 
t.test(result_list.synec.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0007 ********************

# ttest.pct.change.genomes.np.synec.ca09 <- 
t.test(result_list.synec.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.008 ******************************

# ttest.pct.change.genomes.na.synec.ca09 <- 
t.test(result_list.synec.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.015 *********************

# ttest.pct.change.genomes.m.synec.ca09 <- 
t.test(result_list.synec.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.61

# ttest.pct.change.genomes.ns.synec.ca09 <- 
t.test(result_list.synec.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.synec.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.51

#___ca_tolyp___
# normality-negative
# pa
#
# normality-positive
# pb2 pb1 ha np na m ns

# ttest.pct.change.genomes.pb2.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$PB2$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.016 ****************************

# ttest.pct.change.genomes.pb1.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$PB1$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0027 ***********************************

# ttest.pct.change.genomes.ha.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$HA$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.016 *********************************

# ttest.pct.change.genomes.np.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$NP$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.035 ********************************

# ttest.pct.change.genomes.na.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$`NA`$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0011 ********************************

# ttest.pct.change.genomes.m.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$M$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.5

# ttest.pct.change.genomes.ns.tolyp.ca09 <- 
t.test(result_list.tolyp.ca09$NS$pct.change.prop.dvg.sgmt,
       result_list.tolyp.ca09$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.9

#___tx_insu___
# normality-negative
# pb1 ha
#
# normality-positive
# pb2 pa np na m ns

# ttest.pct.change.genomes.pb2.insu.tx12 <- 
t.test(result_list.insu.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.insu.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.15

# ttest.pct.change.genomes.pa.insu.tx12 <- 
t.test(result_list.insu.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.insu.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.22

# ttest.pct.change.genomes.np.insu.tx12 <- 
t.test(result_list.insu.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.insu.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.26

# ttest.pct.change.genomes.na.insu.tx12 <- 
t.test(result_list.insu.tx12$`NA`$pct.change.prop.dvg.sgmt,
       result_list.insu.tx12$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.08

# ttest.pct.change.genomes.m.insu.tx12 <- 
t.test(result_list.insu.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.insu.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.07

# ttest.pct.change.genomes.ns.insu.tx12 <- 
t.test(result_list.insu.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.insu.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.97

#___tx_alpe___
# normality-negative
# NONE
#
# normality-positive
# pb2 pb1 pa ha np na m ns

# ttest.pct.change.genomes.pb2.alpe.tx12 <- 
t.test(result_list.alpe.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.018 *******************************

# ttest.pct.change.genomes.pb1.alpe.tx12 <- 
t.test(result_list.alpe.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.064

# ttest.pct.change.genomes.pa.alpe.tx12 <- 
t.test(result_list.alpe.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.034 ****************************

# ttest.pct.change.genomes.ha.alpe.tx12 <- 
t.test(result_list.alpe.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.64

# ttest.pct.change.genomes.np.alpe.tx12 <- 
t.test(result_list.alpe.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.32

# ttest.pct.change.genomes.na.alpe.tx12 <- 
t.test(result_list.alpe.tx12$`NA`$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.48

# ttest.pct.change.genomes.m.alpe.tx12 <- 
t.test(result_list.alpe.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.14

# ttest.pct.change.genomes.ns.alpe.tx12 <- 
t.test(result_list.alpe.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.alpe.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.46

#___tx_mk2206___
# normality-negative
# pb1 ha np
#
# normality-positive
# pb2 pa na m ns

# ttest.pct.change.genomes.pb2.mk2206.tx12 <- 
t.test(result_list.mk2206.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.mk2206.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.35

# ttest.pct.change.genomes.pa.mk2206.tx12 <- 
t.test(result_list.mk2206.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.mk2206.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.38

# ttest.pct.change.genomes.na.mk2206.tx12 <- 
t.test(result_list.mk2206.tx12$`NA`$pct.change.prop.dvg.sgmt,
       result_list.mk2206.tx12$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.88

# ttest.pct.change.genomes.m.mk2206.tx12 <- 
t.test(result_list.mk2206.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.mk2206.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.013 **************************

# ttest.pct.change.genomes.ns.mk2206.tx12 <- 
t.test(result_list.mk2206.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.mk2206.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.017 *****************************

#___tx_4oi___
# normality-negative
# pb1 
#
# normality-positive
# pb2 pa ha np na m ns

# ttest.pct.change.genomes.pb2.4oi.tx12 <- 
t.test(result_list.4oi.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.32

# ttest.pct.change.genomes.pa.4oi.tx12 <- 
t.test(result_list.4oi.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.27

# ttest.pct.change.genomes.ha.4oi.tx12 <- 
t.test(result_list.4oi.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.237

# ttest.pct.change.genomes.np.4oi.tx12 <- 
t.test(result_list.4oi.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.2325

# ttest.pct.change.genomes.na.4oi.tx12 <- 
t.test(result_list.4oi.tx12$`NA`$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.2372

# ttest.pct.change.genomes.m.4oi.tx12 <- 
t.test(result_list.4oi.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.17

# ttest.pct.change.genomes.ns.4oi.tx12 <- 
t.test(result_list.4oi.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.4oi.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.19

#___tx_uk5099___
# normality-negative
# pb1 
#
# normality-positive
# pb2 pa ha np na m ns

# ttest.pct.change.genomes.pb2.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.086

# ttest.pct.change.genomes.pa.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.34

# ttest.pct.change.genomes.ha.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.21

# ttest.pct.change.genomes.np.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.15

# ttest.pct.change.genomes.pb2.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$`NA`$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.20

# ttest.pct.change.genomes.m.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.15

# ttest.pct.change.genomes.ns.uk5099.tx12 <- 
t.test(result_list.uk5099.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.uk5099.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.21

#___tx_ado___
# normality-negative
# pa na m 
#
# normality-positive
# pb2 pb1 ha np ns

# ttest.pct.change.genomes.pb2.ado.tx12 <- 
t.test(result_list.ado.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.ado.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.009 ******************************

# ttest.pct.change.genomes.pb1.ado.tx12 <- 
t.test(result_list.ado.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.ado.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.002 ******************************

# ttest.pct.change.genomes.ha.ado.tx12 <- 
t.test(result_list.ado.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.ado.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.11

# ttest.pct.change.genomes.np.ado.tx12 <- 
t.test(result_list.ado.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.ado.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.12

# ttest.pct.change.genomes.ns.ado.tx12 <- 
t.test(result_list.ado.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.ado.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.35

#___tx_uri___
# normality-negative
# pa ns ....NOTE: na isn't testable (n=1 DVGs)
#
# normality-positive
# pb2 pb1 ha np m

# ttest.pct.change.genomes.pb2.uri.tx12 <- 
t.test(result_list.uri.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.uri.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.31

# ttest.pct.change.genomes.pb1.uri.tx12 <- 
t.test(result_list.uri.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.uri.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.33

# ttest.pct.change.genomes.ha.uri.tx12 <- 
t.test(result_list.uri.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.uri.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.5

# ttest.pct.change.genomes.np.uri.tx12 <- 
t.test(result_list.uri.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.uri.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.10


# ttest.pct.change.genomes.m.uri.tx12 <- 
t.test(result_list.uri.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.uri.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.31

#___tx_favp___
# normality-negative
# ha na
#
# normality-positive
# pb2 pb1 ha np m

# ttest.pct.change.genomes.pb2.favp.tx12 <- 
t.test(result_list.favp.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.favp.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.11

# ttest.pct.change.genomes.pb1.favp.tx12 <- 
t.test(result_list.favp.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.favp.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.17

# ttest.pct.change.genomes.ha.favp.tx12 <- 
t.test(result_list.favp.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.favp.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.42

# ttest.pct.change.genomes.np.favp.tx12 <- 
t.test(result_list.favp.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.favp.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.09

# ttest.pct.change.genomes.m.favp.tx12 <- 
t.test(result_list.favp.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.favp.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.07

#___tx_lepto___
# normality-negative
# ns.... NOTE: ha + na are not testable(n=2 DVGs)
#
# normality-positive
# pb2 pb1 pa np m

# ttest.pct.change.genomes.pb2.lepto.tx12 <- 
t.test(result_list.lepto.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.lepto.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0003 ****************************

# ttest.pct.change.genomes.pb1.lepto.tx12 <- 
t.test(result_list.lepto.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.lepto.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.080

# ttest.pct.change.genomes.pa.lepto.tx12 <- 
t.test(result_list.lepto.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.lepto.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.016 ********************

# ttest.pct.change.genomes.np.lepto.tx12 <- 
t.test(result_list.lepto.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.lepto.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.15

# ttest.pct.change.genomes.m.lepto.tx12 <- 
t.test(result_list.lepto.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.lepto.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.040 ***********************

#___tx_nosto___
# normality-negative
# ha np na m ns
#
# normality-positive
# pb2 pb1 pa

# ttest.pct.change.genomes.pb2.nosto.tx12 <- 
t.test(result_list.nosto.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.nosto.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.08

# ttest.pct.change.genomes.pb1.nosto.tx12 <- 
t.test(result_list.nosto.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.nosto.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.039 ******************************

# ttest.pct.change.genomes.pa.nosto.tx12 <- 
t.test(result_list.nosto.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.nosto.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.04 *******************************

#___tx_oscil___
# normality-negative
# pb2 na
#
# normality-positive
# pb1 pa ha np m ns

# ttest.pct.change.genomes.pb1.oscil.tx12 <- 
t.test(result_list.oscil.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.oscil.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.19

# ttest.pct.change.genomes.pa.oscil.tx12 <- 
t.test(result_list.oscil.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.oscil.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.06

# ttest.pct.change.genomes.ha.oscil.tx12 <- 
t.test(result_list.oscil.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.oscil.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.11

# ttest.pct.change.genomes.np.oscil.tx12 <- 
t.test(result_list.oscil.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.oscil.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.56

# ttest.pct.change.genomes.m.oscil.tx12 <- 
t.test(result_list.oscil.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.oscil.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.22

# ttest.pct.change.genomes.ns.oscil.tx12 <- 
t.test(result_list.oscil.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.oscil.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.308

#___tx_synec___
# normality-negative
# NONE
#
# normality-positive
# pb2 pb1 pa ha np na m ns

# ttest.pct.change.genomes.pb2.synec.tx12 <- 
t.test(result_list.synec.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.09

# ttest.pct.change.genomes.pb1.synec.tx12 <- 
t.test(result_list.synec.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0059 ***********************

# ttest.pct.change.genomes.pa.synec.tx12 <- 
t.test(result_list.synec.tx12$PA$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$PA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.011 **************************

# ttest.pct.change.genomes.ha.synec.tx12 <- 
t.test(result_list.synec.tx12$HA$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$HA$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.04 *************************

# ttest.pct.change.genomes.np.synec.tx12 <- 
t.test(result_list.synec.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.11

# ttest.pct.change.genomes.na.synec.tx12 <- 
t.test(result_list.synec.tx12$`NA`$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$`NA`$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.06

# ttest.pct.change.genomes.m.synec.tx12 <- 
t.test(result_list.synec.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.083

# ttest.pct.change.genomes.ns.synec.tx12 <- 
t.test(result_list.synec.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.synec.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0849

#___tx_tolyp___
# normality-negative
# pa ha na
#
# normality-positive
# pb2 pb1 np m ns

# ttest.pct.change.genomes.pb2.tolyp.tx12 <- 
t.test(result_list.tolyp.tx12$PB2$pct.change.prop.dvg.sgmt,
       result_list.tolyp.tx12$PB2$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.089

# ttest.pct.change.genomes.pb1.tolyp.tx12 <- 
t.test(result_list.tolyp.tx12$PB1$pct.change.prop.dvg.sgmt,
       result_list.tolyp.tx12$PB1$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.0001 ***************************

# ttest.pct.change.genomes.np.tolyp.tx12 <- 
t.test(result_list.tolyp.tx12$NP$pct.change.prop.dvg.sgmt,
       result_list.tolyp.tx12$NP$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.2

# ttest.pct.change.genomes.m.tolyp.tx12 <- 
t.test(result_list.tolyp.tx12$M$pct.change.prop.dvg.sgmt,
       result_list.tolyp.tx12$M$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.12

# ttest.pct.change.genomes.ns.tolyp.tx12 <- 
t.test(result_list.tolyp.tx12$NS$pct.change.prop.dvg.sgmt,
       result_list.tolyp.tx12$NS$pct.change.count.tvg.sgmt,  
       paired = TRUE, 
       alternative = "two.sided")
#p= 0.3




#  Insu > Alpe > MK2206 > 4-OI > UK5099 > Ado > Uri > Favp > Lepto > Nosto > Oscil > Synec > Tolyp
# We have multiple instances of significance:
# ca09: 
#      insu    np na
#      alpe    -
#      mk2206  m
#      4oi     pa m ns
#      uk5099  ha np na
#      ado     -
#      uri     ha np na
#      favp    pb2 pb1 pa
#      lepto   pb2 pa np na m ns
#      nosto   pb1 pa ha np
#      oscil   pa ha
#      synec   ha np na
#      tolyp   pb2 pb1 ha np na
#      
# tx12:	
#      insu    - 
#      alpe    pb2 pa
#      mk2206  m ns
#      4oi     -
#      uk5099  -
#      ado     pb2 pb1
#      uri     -
#      favp    -
#      lepto   pb2 pa m
#      nosto   pb1 pa
#      oscil   -
#      synec   pb1 pa ha
#      tolyp   pb1
#      
#
#________________________________________________________________________
# //////////////////////////////////  CONTINUE FROM HERE  //////////////////////////////////////
#_______________________________________________________________________________________
# so add significance bars and asterisks
# In this case it's a faceted plot, so I first need to create:
#     an object containing the _*_ for each facet.................  'asterisks_alpevflu1'
#     no need for line segments since the bars of tested groups overlap
#   then my earlier chart object will call this asterisk object.

# asterisk locations
# tx12:	5uM    "PB1", "PA", "HA", "NP", "NA", "M", "NS"
asterisks_alpevflu1 <- data.frame(
  host.ptg = c("MDCK.TX12(H3N2)", "MDCK.TX12(H3N2)", "MDCK.TX12(H3N2)", "MDCK.TX12(H3N2)",
               "MDCK.TX12(H3N2)", "MDCK.TX12(H3N2)", "MDCK.TX12(H3N2)"),              # input facet_grid column names here...include each facet as many times as there are asterisks you need to insert
  geno.sgmt = c("PB1", "PA", "HA", "NP", "NA", "M", "NS"),                               # input facet_grid row names here
  x = c(4, 4, 4, 4, 4, 4, 4),                                                         # X pos'n for each test group with significance
  y = c(50, 50, 50, 50, 50, 50, 50)                                                    # Y pos'n for each test group with significance
)


#__call asterisks into their charts
#
chart.pct.change.tvg_v_prop.dvg.mean.sgmt <-




#______________RESULTS 1.2____________total genomes
# KNOWNS: intimate metabolic cross-talk between flu virus and host metabolism
# UNKNOWN: the impact of alpelisib dosing on flu progeny evolution i.e. SVG vs DVG outcomes
#                   i.   percent change in tvg count ...'count.dvg' 
#                  ii.  percent change in relative abundance of dvgs ... 'prop.dvg'     
#
#
# Ha: %change in prop.dvg ≠ %change in count.tvg
# H0: %percent change in prop.dvg = %change in count.dvg

# Order of events
# i. create percent change object...containing 'pct.change.count.tvg' and 'pct.change.prop.dvg'
# ii. plot chart
# iii. hypothesis-test




#________percent change dataframe
# Okay, here's how you do this
# I want to calculate %change in count.tvg between veh-ctrl(DMSO/dH2O) and the other Tx's 
#      but within strains i.e. for ca09 then tx12
#	'%change' formula = ((new.value - baseline.value) / baseline.value) * 100
#		where
#		- baseline.value   is count.tvg for DMSO/dH2O
#		- new.value	is count.tvg for "Alpe" , "MK2206" , "4-OI" , "UK5099" , "Insu" , "Ado" , "Uri" , "Lepto" , "Nosto" , "Oscil" , "Synec" , "Tolyp"
#
#
#
# To do	new.value - baseline.value ...more accurately count.tvg - count.tvg.baseline
#        I need to create/add a		count.tvg.baseline	variable to the df
#	       but populate it w/the DMSO count.tvg for the DMSO-vehicle grps, and w/dH2O count.tvg for the dH2O-vehicle grps
#	       the challenge is that ea bioreplicate batch has its own unique count.tvg.baseline of DMSO and dH2O
#           so I need to create a df for each bioreplicate(3x) for each strain(2x) and for each vehicle(2x)
#           so     6x dataframes ea for ca09 and tx12...
#		I will do this w/ a combo of    	mutate(...ifelse)     and      fill()            shown in a bit
#
# then bind all 12x df's w/       bind_rows(df1, df2, df3, ... , df6)      6ea for ca09 + tx12
# then plot with ggplot()
# then hypothesis test

# now, make the 12x df's...they'll be bound after in2 the pct.change df
#	here are the df names
# something to keep in mind:
#     these df's are to calculate percent change relative to the vehicle ctrl

#pct.change.genomes.drugvflu1.ca09.dmso.bio1.df
#pct.change.genomes.drugvflu1.ca09.dmso.bio2.df
#pct.change.genomes.drugvflu1.ca09.dmso.bio3.df
#pct.change.genomes.drugvflu1.ca09.dh2o.bio1.df
#pct.change.genomes.drugvflu1.ca09.dh2o.bio2.df
#pct.change.genomes.drugvflu1.ca09.dh2o.bio3.df

#pct.change.genomes.drugvflu1.tx12.dmso.bio1.df
#pct.change.genomes.drugvflu1.tx12.dmso.bio2.df
#pct.change.genomes.drugvflu1.tx12.dmso.bio3.df
#pct.change.genomes.drugvflu1.tx12.dh2o.bio1.df
#pct.change.genomes.drugvflu1.tx12.dh2o.bio2.df
#pct.change.genomes.drugvflu1.tx12.dh2o.bio3.df



### As I create the datafame for ea bioreplicate, 
#       I will assign the unique veh_ctrl(DMSO vs. dH2O) like so:
#
#
# the	   mutate(...ifelse) and fill()	     combo of codes below is to fill-in the new count.tvg.baseline variable for all observations
#	mutate() creates the new variable
#		ifelse() says "if the value in the _Tx_ variable is _DMSO_ i.e. the veh-ctrl, then populate with that observation's count.tvg value
#		if not, populate w/NA
#	then the fill() cmd fills-in the _NA_ values with the last computed value, which the count.tvg of the veh-ctrl (DMSO or dH2O)


library(tidyverse)


#________ca09__dmso__bio1-3__________
pct.change.genomes.drugvflu1.ca09.dmso.bio1.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.CA09(H1N1)" , bioreplicate=="1" , veh.Tx== "DMSO" ) %>%      # we want ca09 dmso bio1
  
  mutate(count.tvg.baseline = ifelse(Tx =="DMSO" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'DMSO' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="DMSO" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="DMSO" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.ca09.dmso.bio2.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.CA09(H1N1)" , bioreplicate=="2" , veh.Tx== "DMSO" ) %>%      # we want ca09 dmso bio2
  
  mutate(count.tvg.baseline = ifelse(Tx =="DMSO" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'DMSO' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="DMSO" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="DMSO" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.ca09.dmso.bio3.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.CA09(H1N1)" , bioreplicate=="3" , veh.Tx== "DMSO" ) %>%      # we want ca09 dmso bio3
  
  mutate(count.tvg.baseline = ifelse(Tx =="DMSO" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'DMSO' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="DMSO" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="DMSO" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

#________ca09__dh2o__bio1-3__________
pct.change.genomes.drugvflu1.ca09.dh2o.bio1.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.CA09(H1N1)" , bioreplicate=="1" , veh.Tx== "dH2O" ) %>%      # we want ca09 dh20 bio1
  
  mutate(count.tvg.baseline = ifelse(Tx =="dH2O" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'dH2O' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="dH2O" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="dH2O" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.ca09.dh2o.bio2.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.CA09(H1N1)" , bioreplicate=="2" , veh.Tx== "dH2O" ) %>%      # we want ca09 dh20 bio2
  
  mutate(count.tvg.baseline = ifelse(Tx =="dH2O" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'dH2O' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="dH2O" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="dH2O" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.ca09.dh2o.bio3.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.CA09(H1N1)" , bioreplicate=="3" , veh.Tx== "dH2O" ) %>%      # we want ca09 dh20 bio3
  
  mutate(count.tvg.baseline = ifelse(Tx =="dH2O" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'dH2O' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="dH2O" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="dH2O" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 





#________tx12__dmso__bio1-3__________
pct.change.genomes.drugvflu1.tx12.dmso.bio1.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.TX12(H3N2)" , bioreplicate=="1" , veh.Tx== "DMSO" ) %>%      # we want tx12 dmso bio1
  
  mutate(count.tvg.baseline = ifelse(Tx =="DMSO" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'DMSO' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="DMSO" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="DMSO" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.tx12.dmso.bio2.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.TX12(H3N2)" , bioreplicate=="2" , veh.Tx== "DMSO" ) %>%      # we want tx12 dmso bio2
  
  mutate(count.tvg.baseline = ifelse(Tx =="DMSO" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'DMSO' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="DMSO" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="DMSO" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.tx12.dmso.bio3.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.TX12(H3N2)" , bioreplicate=="3" , veh.Tx== "DMSO" ) %>%      # we want tx12 dmso bio3
  
  mutate(count.tvg.baseline = ifelse(Tx =="DMSO" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'DMSO' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="DMSO" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="DMSO" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

#________tx12__dh2o__bio1-3__________
pct.change.genomes.drugvflu1.tx12.dh2o.bio1.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.TX12(H3N2)" , bioreplicate=="1" , veh.Tx== "dH2O" ) %>%      # we want tx12 dh20 bio1
  
  mutate(count.tvg.baseline = ifelse(Tx =="dH2O" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'dH2O' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="dH2O" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="dH2O" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.tx12.dh2o.bio2.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.TX12(H3N2)" , bioreplicate=="2" , veh.Tx== "dH2O" ) %>%      # we want tx12 dh20 bio2
  
  mutate(count.tvg.baseline = ifelse(Tx =="dH2O" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'dH2O' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="dH2O" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="dH2O" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 

pct.change.genomes.drugvflu1.tx12.dh2o.bio3.df <- gcay.screen.count %>%       # we're using the object where bioreplicates are not averaged
  filter(host.ptg=="MDCK.TX12(H3N2)" , bioreplicate=="3" , veh.Tx== "dH2O" ) %>%      # we want tx12 dh20 bio3
  
  mutate(count.tvg.baseline = ifelse(Tx =="dH2O" , count.tvg , NA ) ) %>%     # create column for baseline TVG count; allocate the value for vehicle group, but fill the rest of the column with 'NA'
  fill(count.tvg.baseline, .direction = "down") %>%                              # fills-in the 'NA' values with the last computed value, which the count.tvg of the 'dH2O' test group
  mutate(pct.change.count.tvg = ((count.tvg - count.tvg.baseline)/count.tvg.baseline) * 100) %>%    # create a variable to hold the percent change of TVG count between the vehicle-treated group and each drug/metabolite/biologic treatment 
  mutate(pct.change.count.tvg = round(pct.change.count.tvg, digits =2)       #I'm rounding to two decimal places because I'll be pasting these values on the plotted chart
  ) %>%
  
  mutate(prop.dvg.baseline = ifelse(Tx =="dH2O" , prop.dvg , NA ) ) %>%     # same as count.tvg above, but for prop.dvg
  fill(prop.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.prop.dvg = ((prop.dvg - prop.dvg.baseline)/prop.dvg.baseline) * 100) %>%  
  mutate(pct.change.prop.dvg = round(pct.change.prop.dvg, digits =2)      
  ) %>%
  
  mutate(count.dvg.baseline = ifelse(Tx =="dH2O" , count.dvg , NA ) ) %>%     # same as count.tvg above, but for count.dvg
  fill(count.dvg.baseline, .direction = "down") %>%                              
  mutate(pct.change.count.dvg = ((count.dvg - count.dvg.baseline)/count.dvg.baseline) * 100) %>%  
  mutate(pct.change.count.dvg = round(pct.change.count.dvg, digits =2)      
  ) 



#_____bind %change data_____
# here I bind the earlier pct.change.pakt.au.mock data from all bioreplicates into a single df
# I will then filter the records as needed to conduct anovas or plot charts

pct.change.genomes.drugvflu1.bind <- bind_rows(pct.change.genomes.drugvflu1.ca09.dmso.bio1.df,
                                              pct.change.genomes.drugvflu1.ca09.dmso.bio2.df,
                                              pct.change.genomes.drugvflu1.ca09.dmso.bio3.df,
                                              pct.change.genomes.drugvflu1.ca09.dh2o.bio1.df,
                                              pct.change.genomes.drugvflu1.ca09.dh2o.bio2.df,
                                              pct.change.genomes.drugvflu1.ca09.dh2o.bio3.df,
                                              pct.change.genomes.drugvflu1.tx12.dmso.bio1.df,
                                              pct.change.genomes.drugvflu1.tx12.dmso.bio2.df,
                                              pct.change.genomes.drugvflu1.tx12.dmso.bio3.df,
                                              pct.change.genomes.drugvflu1.tx12.dh2o.bio1.df,
                                              pct.change.genomes.drugvflu1.tx12.dh2o.bio2.df,
                                              pct.change.genomes.drugvflu1.tx12.dh2o.bio3.df)

#View(pct.change.genomes.drugvflu1.bind)

#____________________barchart________________________________
# Here, I plot charts 
#     both point charts and barcharts

# first we'll show:
# pct.change.count.tvg and pct.change.prop.dvg...full spread of vehicle-corrected bioreplicates
chart.pct.change.tvg_v_prop.dvg <-
  ggplot(pct.change.genomes.drugvflu1.bind, aes(x = Tx, shape = bioreplicate)) +
  geom_point(aes(y = pct.change.count.tvg, color = "pct.change.count.tvg", shape= bioreplicate), position = position_jitter(width = 0.05), size = 3, alpha = 0.5) +
  geom_point(aes(y = pct.change.prop.dvg, color = "pct.change.prop.dvg", shape= bioreplicate), position = position_jitter(width = 0.05), size = 3, alpha = 0.5) +
  labs(y = "% Change", x= "Treatment") +
  #geom_text(aes(y = pct.change.count.tvg, label = round(pct.change.count.tvg, 1)), vjust = -0.5, size = 2.5, nudge_x = -0.2) +     # an alternate to nudge_x is      position = position_jitter(width = 0.4)
  #geom_text(aes(y = pct.change.prop.dvg, label = round(pct.change.prop.dvg, 1)), vjust = -0.5, size = 2.5, nudge_x = 0.2) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(title = "Bioreplicate", direction = "horizontal"),
         color = guide_legend(title = "Variable", direction = "horizontal")) +
  scale_color_manual(values = c("pct.change.count.tvg" = "blue", "pct.change.prop.dvg" = "red")) +
  scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 18)) +
  facet_wrap(~host.ptg,
             nrow = 2,
             scales = "free_y")


# Next we'll show:
# mean.pct.change.tcu.titr and mean.pct.change.prop.ncu ... i.e. bioreplicates aggregated
chart.pct.change.tvg_v_prop.dvg.mean <- 
  pct.change.genomes.drugvflu1.bind %>%
  group_by(host.ptg, Tx, veh.Tx       # grouping for mean of bioreplicates
  ) %>% 
  summarise(n = n(),
            mean.pct.tvg = mean(pct.change.count.tvg),     # mean TVG count (dosing experiment bioreplicates)
            sd.tvg = sd(pct.change.count.tvg),               # sd
            sem.tvg = sd.tvg /sqrt(n()),          #sem
            
            mean.pct.dvg = mean(pct.change.count.dvg),     # mean DVG count (dosing experiment bioreplicates)
            sd.dvg = sd(pct.change.count.dvg),               #sd
            sem.dvg = sd.dvg /sqrt(n()),          #sem
            
            mean.prop.dvg = mean(pct.change.prop.dvg),     # mean DVG relative abundance (dosing experiment bioreplicates)
            sd.prop.dvg = sd(pct.change.prop.dvg),           #sd
            sem.prop.dvg = sd.prop.dvg /sqrt(n())   #sem
  ) %>%
  ungroup() %>%
  rename(
    pct.change.count.tvg = mean.pct.tvg,                                           # rename variable 
    pct.change.count.dvg = mean.pct.dvg,
    pct.change.prop.dvg = mean.prop.dvg,                                           # rename variable
  ) %>% 
  ggplot(aes(x = Tx)) +
  geom_bar(aes(y = pct.change.count.tvg, fill = "pct.change.count.tvg"), position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), stat = "identity", alpha = 0.5) +
  geom_errorbar(aes(y = pct.change.count.tvg, ymin = pct.change.count.tvg - sem.tvg, ymax = pct.change.count.tvg),    # some error bars are too long, so ymax is configured to not show
                position = position_dodge(width = 0.9), width = 0.2, alpha = 0.5) +
  geom_bar(aes(y = pct.change.prop.dvg, fill = "pct.change.prop.dvg"), position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), stat = "identity", alpha = 0.5) +
  geom_errorbar(aes(y = pct.change.prop.dvg, ymin = pct.change.prop.dvg - sem.prop.dvg, ymax = pct.change.prop.dvg + sem.prop.dvg),   # error bars are manageable, so the full range is shown
                position = position_dodge(width = 0.9), width = 0.2, alpha = 0.5) +
  labs(y = "% Change", x= "Treatment", fill = NULL) +  # Remove legend title with 'fill=NULL'...it's clunky and the 
  geom_text(aes(y = pct.change.count.tvg, label = round(pct.change.count.tvg, 1)), angle=0,  vjust = -0.5, size = 2.5, nudge_x = -0.2) +     # an alternate to nudge_x is      position = position_jitter(width = 0.4)
  geom_text(aes(y = pct.change.prop.dvg, label = round(pct.change.prop.dvg, 1)), angle=0, vjust = -0.5, size = 2.5, nudge_x = 0.2) +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=7)) +
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("pct.change.count.tvg" = "blue", "pct.change.prop.dvg" = "red")) +
  facet_wrap(~host.ptg,
             nrow = 2,
             scales = "free_y")  

# NEXT do the hypothesis testing to determine significance,
#    then add the asterisks and line segments denoting significance to the chart object





#____________Hypothesis Testing_________________
# we are testing percent change of count.tvg against percent change of prop.dvg within each sample (paired)
# percent change is normalized against vehicle-treatment control so averages can be taken with relative safety (parametric t-test), 
#      unless data set is normality-negative(non-parametric rank test)
#
#
# for each variable you want to test, create and object of its values from all 3x bioreplicates
#

#___promising Tx's
#  ca09:    4oi uk5099 alpe lepto
#  tx12:    4oi uk5099 favp ado lepto
#_____ca09___pct.change.count.tvg
pct.change.genomes.tvg.4oi.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="4-OI") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.uk5099.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="UK5099") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.alpe.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="Alpe") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.lepto.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="Lepto") %>%
  select(pct.change.count.tvg)



#_____ca09___pct.change.prop.dvg
pct.change.genomes.prop.dvg.4oi.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="4-OI") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.uk5099.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="UK5099") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.alpe.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="Alpe") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.lepto.ca09 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.CA09(H1N1)" , Tx =="Lepto") %>%
  select(pct.change.prop.dvg)



#  tx12:    4oi uk5099 favp ado lepto
#_____tx12___pct.change.count.tvg
pct.change.genomes.tvg.4oi.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="4-OI") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.uk5099.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="UK5099") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.favp.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="Favp") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.ado.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="Ado") %>%
  select(pct.change.count.tvg)
pct.change.genomes.tvg.lepto.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="Lepto") %>%
  select(pct.change.count.tvg)


#_____tx12___pct.change.prop.dvg
pct.change.genomes.prop.dvg.4oi.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="4-OI") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.uk5099.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="UK5099") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.favp.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="Favp") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.ado.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="Ado") %>%
  select(pct.change.prop.dvg)
pct.change.genomes.prop.dvg.lepto.tx12 <- pct.change.genomes.drugvflu1.bind %>%
  filter( host.ptg=="MDCK.TX12(H3N2)" , Tx =="Lepto") %>%
  select(pct.change.prop.dvg)






#___________normality
#test normality one group at a time
# groups with promising separation are:
#  ca09:    4oi uk5099 alpe lepto
#  tx12:    4oi uk5099 favp ado lepto

# ca09
shapiro.test(pct.change.genomes.prop.dvg.4oi.ca09$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.4oi.ca09$pct.change.count.tvg)
#normality-positive  (p > 0.05)

shapiro.test(pct.change.genomes.prop.dvg.uk5099.ca09$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.uk5099.ca09$pct.change.count.tvg)
#normality-positive  (p > 0.05)

shapiro.test(pct.change.genomes.prop.dvg.alpe.ca09$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.alpe.ca09$pct.change.count.tvg)
#normality-positive  (p > 0.05)

shapiro.test(pct.change.genomes.prop.dvg.lepto.ca09$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.lepto.ca09$pct.change.count.tvg)
#normality-positive  (p > 0.05)



#  tx12:    4oi uk5099 favp ado lepto
# tx12
shapiro.test(pct.change.genomes.prop.dvg.4oi.tx12$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.4oi.tx12$pct.change.count.tvg)
#normality-positive  (p > 0.05)

shapiro.test(pct.change.genomes.prop.dvg.uk5099.tx12$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.uk5099.tx12$pct.change.count.tvg)
#normality-positive  (p > 0.05)

# tx12_favp
shapiro.test(pct.change.genomes.prop.dvg.favp.tx12$pct.change.prop.dvg)  # normality-negative (p < 0.05)
shapiro.test(pct.change.genomes.tvg.favp.tx12$pct.change.count.tvg)
#normality-positive  (p > 0.05)

shapiro.test(pct.change.genomes.prop.dvg.ado.tx12$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.ado.tx12$pct.change.count.tvg)
#normality-positive  (p > 0.05)

shapiro.test(pct.change.genomes.prop.dvg.lepto.tx12$pct.change.prop.dvg)
shapiro.test(pct.change.genomes.tvg.lepto.tx12$pct.change.count.tvg)
#normality-positive  (p > 0.05)





#__________equality of variances
# test equality of variance
# Not necessary here:
#   equal variances are not a required assumption for paired T-Test, only that data is normality-positive



#_____parametric or non-parametric test?
#
# normality-negative data sets will undergo Wilcoxon Signed Rank Test
# normality-positive, variance-unequal data sets will undergo Welch's T-Test
# normality-positive, variance-equal data sets will undergo Student's T-Test
# ALL tests are paired, since test variables are from the same sample

#_______Paired one-sample ttest

# ca09:    4oi uk5099 alpe lepto
ttest.pct.change.genomes.4oi.ca09 <- 
  t.test(pct.change.genomes.prop.dvg.4oi.ca09$pct.change.prop.dvg, 
         pct.change.genomes.tvg.4oi.ca09$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.003 **********

ttest.pct.change.genomes.uk5099.ca09 <- 
  t.test(pct.change.genomes.prop.dvg.uk5099.ca09$pct.change.prop.dvg, 
         pct.change.genomes.tvg.uk5099.ca09$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.056  .... Really God?!

ttest.pct.change.genomes.alpe.ca09 <- 
  t.test(pct.change.genomes.prop.dvg.alpe.ca09$pct.change.prop.dvg, 
         pct.change.genomes.tvg.alpe.ca09$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.18

ttest.pct.change.genomes.lepto.ca09 <- 
  t.test(pct.change.genomes.prop.dvg.lepto.ca09$pct.change.prop.dvg, 
         pct.change.genomes.tvg.lepto.ca09$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.13


# tx12:    4oi uk5099 favp ado lepto

#_______Paired one-sample ttest
ttest.pct.change.genomes.4oi.tx12 <- 
  t.test(pct.change.genomes.prop.dvg.4oi.tx12$pct.change.prop.dvg, 
         pct.change.genomes.tvg.4oi.tx12$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = FALSE,        
         alternative = "two.sided")
# p=0.25

ttest.pct.change.genomes.uk5099.tx12 <- 
  t.test(pct.change.genomes.prop.dvg.uk5099.tx12$pct.change.prop.dvg, 
         pct.change.genomes.tvg.uk5099.tx12$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = FALSE,        
         alternative = "two.sided")
# p=0.23

ttest.pct.change.genomes.favp.tx12 <- 
  t.test(pct.change.genomes.prop.dvg.favp.tx12$pct.change.prop.dvg, 
         pct.change.genomes.tvg.favp.tx12$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.037 ***********

ttest.pct.change.genomes.ado.tx12 <- 
  t.test(pct.change.genomes.prop.dvg.ado.tx12$pct.change.prop.dvg, 
         pct.change.genomes.tvg.ado.tx12$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.01 ************

ttest.pct.change.genomes.lepto.tx12 <- 
  t.test(pct.change.genomes.prop.dvg.lepto.tx12$pct.change.prop.dvg, 
         pct.change.genomes.tvg.lepto.tx12$pct.change.count.tvg,
         paired = TRUE,               # One sample, paired observations
         var.equal = TRUE,        
         alternative = "two.sided")
# p=0.13


# We have multiple instances of significance:
#    ca09: 4oi
#    tx12: ado, favp
#
# so add significance bars and asterisks
# In this case it's a faceted plot, so I first need to create:
#     an object containing the _*_ for each facet.................  'asterisks_drugvflu1'
#     no need for line segments since the bars of tested groups overlap
#   then my earlier chart object will call this asterisk object.


asterisks_drugvflu1 <- data.frame(
  host.ptg = c("MDCK.CA09(H1N1)", "MDCK.TX12(H3N2)", "MDCK.TX12(H3N2)" ),  # input facet values here...include each facet as many times as there are asterisks you need to insert
  x = c(6, 8, 10),                                                         # X pos'n for each facet
  y = c(375, 90, 250)                                                    # Y pos'n for each facet
)

chart.pct.change.tvg_v_prop.dvg.mean <- 
  chart.pct.change.tvg_v_prop.dvg.mean +
  geom_text(data = asterisks_drugvflu1,
            aes(x = x, y = y), label = "*", size = 5)                 #sig _*_ ...calls data from 




#_______export the plot...conveniently export plots w/ggsave()
# Set desired save directory
save_directory <- "~/R_gcay_screen"
# save to above directory
ggsave(file.path(save_directory, "chart.pct.change.tvg_v_prop.dvg.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       chart.pct.change.tvg_v_prop.dvg, 
       width = 7 , height = 9)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements

ggsave(file.path(save_directory, "chart.pct.change.tvg_v_prop.dvg.mean.tiff"),    # "" contains desired filename...MUST include file extension...here its .tiff
       chart.pct.change.tvg_v_prop.dvg.mean, 
       width = 7 , height = 9)         # width+height numbers are in inches...6x6 works for a standard 8x11 sheet...ultimately depends on journal requirements


