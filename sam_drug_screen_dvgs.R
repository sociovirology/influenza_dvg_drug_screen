library(ggplot2)
library(forcats)
library(dplyr)
library(tidyr)
library(multcomp)
library(lme4)
library(nlme)
library(ggthemes)
library(emmeans)
library(ggbeeswarm)


#This code below assumes I've run 0.1, 0.2 and result1.1.

#Overall
head(gcay.screen.count)
#Averages
head(gcay_screen_bio.avg)

#Per Segment
head(gcay.screen.count.sgmt)
#Per segment averages
head(gcay_screen_bio.avg.sgmt)


### Core results
#First, let's do counts of DVGs. And let's first examine data

#Counts
ggplot(gcay.screen.count, aes(x=Tx, y=count.dvg, color=host.ptg)) +
  geom_point() + theme(legend.position = "none") + facet_wrap(~host.ptg) 

#Are DVG counts correlated to total counts?
ggplot(gcay.screen.count, aes(x=count.tvg, y=count.dvg, fill=host.ptg)) +
  geom_point() + theme(legend.position = "none") + facet_wrap(~host.ptg)
#Pretty related, but divergence in the top range

#Summary stats
#For all treatments
subset(gcay_screen_bio.avg, select = c("Tx","mean.count.tvg", "mean.count.dvg", "sd.dvg", "mean.prop.dvg"))

#Proportions
ggplot(gcay.screen.count, aes(x=Tx, y=prop.dvg, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(bioreplicate~host.ptg)

ggplot(gcay.screen.count, aes(x=Tx, y=prop.dvg, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(~host.ptg)

#Facet by Vehicle DVG
ggplot(gcay.screen.count, aes(x=Tx, y=prop.dvg, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~veh.Tx, scales = "free_x")

#Facet by Vehicle TVG
ggplot(gcay.screen.count, aes(x=Tx, y=count.tvg, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~veh.Tx, scales = "free_x")


###### Stat Testing ######
#Seems pretty straightforward to do an ANOVA

#First curious if any of them have a consistent effect regardless of strain
model0_all <- lm(prop.dvg ~ Tx, data = gcay.screen.count)
summary(model0_all)
summary.aov(model0_all)
#NS

#Control for bioreplicate
model1_all <- lm(prop.dvg ~ Tx + bioreplicate, data = gcay.screen.count)
summary(model1_all)
summary.aov(model1_all)
#NS

#Control for bioreplicate and strain 
model2_all <- lm(prop.dvg ~ Tx + bioreplicate + host.ptg, data = gcay.screen.count)
summary(model2_all)
summary.aov(model2_all)
#Host very different and significant, so split

#Wait! Control for vehicle
model3_all <- lm(prop.dvg ~ Tx + veh.Tx, data = gcay.screen.count)
summary(model3_all)
summary.aov(model3_all)
#NS

### CA09 First
## Prop DVG
model0_ca <- lm(prop.dvg ~ Tx, data = subset(gcay.screen.count, host.ptg == "MDCK.CA09(H1N1)"))
summary(model0_ca)
summary.aov(model0_ca)
#NS

#Control for bioreplicate
model1_ca <- lm(prop.dvg ~ Tx + bioreplicate, data = subset(gcay.screen.count, host.ptg == "MDCK.CA09(H1N1)"))
summary(model1_ca)
summary.aov(model1_ca)
#NS

#Control for bioreplicate and vehicle
model2_ca <- lm(prop.dvg ~ Tx + bioreplicate + veh.Tx, data = subset(gcay.screen.count, host.ptg == "MDCK.CA09(H1N1)"))
summary(model2_ca)
summary.aov(model2_ca)
#NS

#Let's split up data frame into host and vehicle
dvg_overall_ca_dh20 <-  gcay.screen.count %>%
  filter(host.ptg == "MDCK.CA09(H1N1)") %>%
  filter(veh.Tx == "dH2O")

dvg_overall_ca_dmso <-  gcay.screen.count %>%
  filter(host.ptg == "MDCK.CA09(H1N1)") %>%
  filter(veh.Tx == "DMSO")

### CA09 dH20
#Prop DVG
model0_dvg_overall_ca_dh20 <- lm(prop.dvg ~ Tx, data = dvg_overall_ca_dh20)
summary(model0_dvg_overall_ca_dh20)
summary.aov(model0_dvg_overall_ca_dh20)
#NS

#Prop DVG control for bioreplicate
model1_dvg_overall_ca_dh20 <- lm(prop.dvg ~ Tx + bioreplicate, data = dvg_overall_ca_dh20)
summary(model1_dvg_overall_ca_dh20)
summary.aov(model1_dvg_overall_ca_dh20)
#Overall borderline, bioreplicate significant

#Proportion of Variance Explained
(summary.aov(model1_dvg_overall_ca_dh20)[[1]]$'Sum Sq' / sum(summary.aov(model1_dvg_overall_ca_dh20)[[1]]$'Sum Sq'))*100
#     Tx          bioreplicate  Residuals 
#[1] 21.63242 35.55613 42.81145


#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(model1_dvg_overall_ca_dh20, las = 1)
par(opar)
#very nice

#Multiple comparisons with a control model 1 
model1_dvg_overall_ca_dh20.glht.ctrl <- glht(model1_dvg_overall_ca_dh20, linfct = mcp(Tx= "Dunnet"))
summary(model1_dvg_overall_ca_dh20.glht.ctrl)
#Linear Hypotheses:
#  Estimate Std. Error t value Pr(>|t|)
#Insu - dH2O == 0  -1.377e-04  1.351e-03  -0.102    1.000
#Ado - dH2O == 0   -8.786e-04  1.351e-03  -0.650    0.985
#Uri - dH2O == 0   -7.352e-04  1.351e-03  -0.544    0.995
#Lepto - dH2O == 0  5.197e-05  1.351e-03   0.038    1.000
#Nosto - dH2O == 0  4.936e-04  1.351e-03   0.365    1.000
#Oscil - dH2O == 0 -2.251e-03  1.351e-03  -1.666    0.457
#Synec - dH2O == 0 -1.827e-04  1.351e-03  -0.135    1.000
#Tolyp - dH2O == 0 -2.083e-03  1.351e-03  -1.542    0.535
#(Adjusted p values reported -- single-step method)

#Nope!

#Now, simultaneous confidence intervals
confint(model1_dvg_overall_ca_dh20.glht.ctrl)
plot(confint(model1_dvg_overall_ca_dh20.glht.ctrl))


### CA09 DMSO
#Prop DVG
model0_dvg_overall_ca_dmso <- lm(prop.dvg ~ Tx, data = dvg_overall_ca_dmso)
summary(model0_dvg_overall_ca_dmso)
summary.aov(model0_dvg_overall_ca_dmso)
#NS

#Prop DVG control for bioreplicate
model1_dvg_overall_ca_dmso <- lm(prop.dvg ~ Tx + bioreplicate, data = dvg_overall_ca_dmso)
summary(model1_dvg_overall_ca_dmso)
summary.aov(model1_dvg_overall_ca_dmso)
#Still NS

#Multiple comparisons with a control model 1 
model1_dvg_overall_ca_dmso.glht.ctrl <- glht(model1_dvg_overall_ca_dmso, linfct = mcp(Tx= "Dunnet"))
summary(model1_dvg_overall_ca_dmso.glht.ctrl)
#Linear Hypotheses:
#  Estimate Std. Error t value Pr(>|t|)
#Alpe - DMSO == 0   -0.0018475  0.0025124  -0.735    0.912
#MK2206 - DMSO == 0 -0.0033345  0.0025124  -1.327    0.568
#4-OI - DMSO == 0   -0.0014392  0.0025124  -0.573    0.965
#UK5099 - DMSO == 0 -0.0002803  0.0025124  -0.112    1.000
#Favp - DMSO == 0    0.0003818  0.0025124   0.152    1.000

#Nothing

#Now, simultaneous confidence intervals
confint(model1_dvg_overall_ca_dmso.glht.ctrl)
plot(confint(model1_dvg_overall_ca_dmso.glht.ctrl))
#Nope!

### TVG ALL
model0_TVG_ca <- lm(count.tvg ~ Tx, data = subset(gcay.screen.count, host.ptg == "MDCK.CA09(H1N1)"))
summary(model0_TVG_ca)
summary.aov(model0_TVG_ca)
#NS

#Control for bioreplicate
model1_TVG_ca <- lm(count.tvg ~ Tx + bioreplicate, data = subset(gcay.screen.count, host.ptg == "MDCK.CA09(H1N1)"))
summary(model1_TVG_ca)
summary.aov(model1_TVG_ca)
#NS

### TVG by Vehicle dH20
model0_TVG_overall_ca_dh20 <- lm(count.tvg ~ Tx, data = dvg_overall_ca_dh20)
summary(model0_TVG_overall_ca_dh20)
summary.aov(model0_TVG_overall_ca_dh20)
#NS

# TVG by Vehicle dH20 with bioreplicate control
model1_TVG_overall_ca_dh20 <- lm(count.tvg ~ Tx + bioreplicate, data = dvg_overall_ca_dh20)
summary(model1_TVG_overall_ca_dh20)
summary.aov(model1_TVG_overall_ca_dh20)
#NS

### TVG by Vehicle DMSO
model0_TVG_overall_ca_dmso <- lm(count.tvg ~ Tx, data = dvg_overall_ca_dmso)
summary(model0_TVG_overall_ca_dmso)
summary.aov(model0_TVG_overall_ca_dmso)
#NS

# TVG by Vehicle DMSO with bioreplicate control
model1_TVG_overall_ca_dmso <- lm(count.tvg ~ Tx + bioreplicate, data = dvg_overall_ca_dmso)
summary(model1_TVG_overall_ca_dmso)
summary.aov(model1_TVG_overall_ca_dmso)
#Significant, hit with 4-OI 

#Proportion of variance explained
(summary.aov(model1_TVG_overall_ca_dmso)[[1]]$'Sum Sq' / sum(summary.aov(model1_TVG_overall_ca_dmso)[[1]]$'Sum Sq'))*100
# Tx        bioreplicate  Residuals
# 40.33672  29.15572      30.50756
#Nice chunk of variance explained! 

#Check Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(model1_TVG_overall_ca_dmso, las = 1)
par(opar)
#Okay data point 18 is a bit off

#Multiple comparisons with a control model 1 
model1_TVG_overall_ca_dmso.glht.ctrl <- glht(model1_TVG_overall_ca_dmso, linfct = mcp(Tx= "Dunnet"))
summary(model1_TVG_overall_ca_dmso.glht.ctrl)
#Linear Hypotheses:
#Estimate Std. Error t value Pr(>|t|)  
#Alpe - DMSO == 0      32675      38086   0.858   0.8550  
#MK2206 - DMSO == 0   -18943      38086  -0.497   0.9802  
#4-OI - DMSO == 0     103797      38086   2.725   0.0769 .
#UK5099 - DMSO == 0    44272      38086   1.162   0.6732  
#Favp - DMSO == 0      60615      38086   1.592   0.4149 

#4-OI now drops out as borderline

#Now, simultaneous confidence intervals
confint(model1_TVG_overall_ca_dmso.glht.ctrl)
plot(confint(model1_TVG_overall_ca_dmso.glht.ctrl))
#Nope!

### CA09 End

### TX12 Second
## Prop DVG
model0_tx <- lm(prop.dvg ~ Tx, data = subset(gcay.screen.count, host.ptg == "MDCK.TX12(H3N2)"))
summary(model0_tx)
summary.aov(model0_tx)
#Significant right out of the gate!
#Hits on Insulin, MK2206, Ado, Uri

#Control for bioreplicate
model1_tx <- lm(prop.dvg ~ Tx + bioreplicate, data = subset(gcay.screen.count, host.ptg == "MDCK.TX12(H3N2)"))
summary(model1_tx)
summary.aov(model1_tx)
#Still significant, same ones

#Let's split up data frame into host and vehicle
dvg_overall_tx_dh20 <-  gcay.screen.count %>%
  filter(host.ptg == "MDCK.TX12(H3N2)") %>%
  filter(veh.Tx == "dH2O")

dvg_overall_tx_dmso <-  gcay.screen.count %>%
  filter(host.ptg == "MDCK.TX12(H3N2)") %>%
  filter(veh.Tx == "DMSO")

### TX12 dH20
#Prop DVG
model0_dvg_overall_tx_dh20 <- lm(prop.dvg ~ Tx, data = dvg_overall_tx_dh20)
summary(model0_dvg_overall_tx_dh20)
summary.aov(model0_dvg_overall_tx_dh20)
#Significant hits on Insulin, Ado, Uri
#borderline on Lepto  Nosto Oscil

#Prop DVG control for bioreplicate
model1_dvg_overall_tx_dh20 <- lm(prop.dvg ~ Tx + bioreplicate, data = dvg_overall_tx_dh20)
summary(model1_dvg_overall_tx_dh20)
summary.aov(model1_dvg_overall_tx_dh20)
#Same picture

#Proportion of Variance Explained
(summary.aov(model1_dvg_overall_tx_dh20)[[1]]$'Sum Sq' / sum(summary.aov(model1_dvg_overall_tx_dh20)[[1]]$'Sum Sq'))*100
#     Tx          bioreplicate  Residuals 
#[1]  65.980860  4.105028      29.914111
#65% of variance, that's so high!

#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(model1_dvg_overall_tx_dh20, las = 1)
par(opar)
#beautiful!

#Multiple comparisons with a control model 1 
model1_dvg_overall_tx_dh20.glht.ctrl <- glht(model1_dvg_overall_tx_dh20, linfct = mcp(Tx= "Dunnet"))
summary(model1_dvg_overall_tx_dh20.glht.ctrl)
#Insu - dH2O == 0  0.006694   0.002327   2.877   0.0603 .  
#Ado - dH2O == 0   0.011644   0.002327   5.004   <0.001 ***
#  Uri - dH2O == 0   0.007205   0.002327   3.096   0.0395 *  
#  Lepto - dH2O == 0 0.004718   0.002327   2.027   0.2698    
#Nosto - dH2O == 0 0.004279   0.002327   1.839   0.3589    
#Oscil - dH2O == 0 0.004587   0.002327   1.971   0.2945    
#Synec - dH2O == 0 0.001646   0.002327   0.708   0.9763    
#Tolyp - dH2O == 0 0.002348   0.002327   1.009   0.8727

#Ado and Uri there, Insulin drops off

#Now, simultaneous confidence intervals
confint(model1_dvg_overall_tx_dh20.glht.ctrl)
plot(confint(model1_dvg_overall_tx_dh20.glht.ctrl))



### TX12 DMSO
#Prop DVG
model0_dvg_overall_tx_dmso <- lm(prop.dvg ~ Tx, data = dvg_overall_tx_dmso)
summary(model0_dvg_overall_tx_dmso)
summary.aov(model0_dvg_overall_tx_dmso)
#Now NS overall, but still hit on MK2206

#Prop DVG control for bioreplicate
model1_dvg_overall_tx_dmso <- lm(prop.dvg ~ Tx + bioreplicate, data = dvg_overall_tx_dmso)
summary(model1_dvg_overall_tx_dmso)
summary.aov(model1_dvg_overall_tx_dmso)
#Still NS, MK2206 now borderline 

#Multiple comparisons with a control model 1 
model1_dvg_overall_tx_dmso.glht.ctrl <- glht(model1_dvg_overall_tx_dmso, linfct = mcp(Tx= "Dunnet"))
summary(model1_dvg_overall_tx_dmso.glht.ctrl)
#Linear Hypotheses:
#Estimate Std. Error t value Pr(>|t|)
#Alpe - DMSO == 0   -3.966e-04  3.621e-03  -0.110    1.000
#MK2206 - DMSO == 0  7.909e-03  3.621e-03   2.185    0.179
#4-OI - DMSO == 0    7.838e-04  3.621e-03   0.216    1.000
#UK5099 - DMSO == 0 -7.570e-05  3.621e-03  -0.021    1.000
#Favp - DMSO == 0    2.963e-05  3.621e-03   0.008    1.000

#Nada

#Now, simultaneous confidence intervals
confint(model1_dvg_overall_tx_dmso.glht.ctrl)
plot(confint(model1_dvg_overall_tx_dmso.glht.ctrl))
#Nada

## TVG
model0_TVG_tx <- lm(count.tvg ~ Tx, data = subset(gcay.screen.count, host.ptg == "MDCK.TX12(H3N2)"))
summary(model0_TVG_tx)
summary.aov(model0_TVG_tx)
#NS

#Control for bioreplicate
model1_TVG_tx <- lm(count.tvg ~ Tx + bioreplicate, data = subset(gcay.screen.count, host.ptg == "MDCK.TX12(H3N2)"))
summary(model1_TVG_tx)
summary.aov(model1_TVG_tx)
#NS

### TX12 End

## Per Segment Stats ##
#First let's overview the data
#Proportions
ggplot(gcay.screen.count.sgmt, aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx)
#Informative
#Looking at some of those proportions, it seems obvious that the count of segments is super important

#Segment Count
ggplot(gcay.screen.count.sgmt, aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx)
#Seems like a lot of increases in non-polymerase segments, but hard to appreciate polymerase segments

#Polymerase segments only
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#Other segments
non_polymerase_segments <- c("HA", "NP", "NA", "M", "NS")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  non_polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#Now let's visualize proportions the same way
#Prop DVG Polymerase segments only
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#Prop DVG Other segments
non_polymerase_segments <- c("HA", "NP", "NA", "M", "NS")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  non_polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#Facet by Vehicle DVG
ggplot(gcay.screen.count, aes(x=Tx, y=prop.dvg, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~veh.Tx, scales = "free_x")

#Facet by Vehicle TVG
ggplot(gcay.screen.count, aes(x=Tx, y=count.tvg, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~veh.Tx, scales = "free_x")

#Let's start with simple models that are probably not even very useful
model0_segment_dvg_all <- lm(prop.dvg.sgmt ~ Tx, data = gcay.screen.count.sgmt)
summary(model0_segment_dvg_all)
summary.aov(model0_segment_dvg_all)
#Multiple R-squared:  0.1312,	Adjusted R-squared:  0.1139 
#F-statistic: 7.571 on 14 and 702 DF,  p-value: 5.792e-15
#Significant, hits on Insu, MK2206, Ado and Uri. Alpe borderline

#Just to have a visual, this is what this looks like
ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(~host.ptg)
#This plot is actually super, super informative!!!!

#Another view not accounting for strain
ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate, shape=host.ptg)) +
  geom_point() + theme(legend.position = "none")
#Publication Settings
figure2a <- ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  xlab("Treatment") +
  ylab("Proportion of Defective Viral Genomes (DelVGs)") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 14)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  scale_shape_manual(values=c(21, 23)) +
  geom_point(aes(fill = bioreplicate, shape = host.ptg), size = 3, colour = "black", position=position_jitter(width=0.2), alpha=0.90) +
  scale_fill_colorblind() +
  #stat_summary(geom="errorbar", fun = "mean", color="black", size=5) +
  #stat_summary(aes(yintercept = after_stat(y), x = 0), fun = mean, geom = "hline", orientation = "x") +
  #geom_hline(yintercept = mean(gcay.screen.count.sgmt$prop.dvg.sgmt)) +
  #geom_hline(yintercept = mean) +
  #geom_errorbar(size=0.8, width=0.2, aes(ymax=..y..,ymin=..y.., y = mean)) +
  theme(legend.position = "none")
#Diamonds are TX12, circles CA09
figure2a
ggsave(figure2a, filename = "figure2a.pdf")

#Another view not accounting for strain now Total Viral Genomes
ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate, shape=host.ptg)) +
  geom_point() + theme(legend.position = "none")

figure2b <- ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate)) +
  xlab("Treatment") +
  ylab("Total Viral Genomes") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 14)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  #geom_quasirandom() +
  scale_shape_manual(values=c(21, 23)) +
  geom_point(aes(fill = bioreplicate, shape = host.ptg), size = 3, colour = "black", position=position_jitter(width=0.2), alpha=0.90) +
  scale_fill_colorblind() +
  #geom_hline(yintercept = mean) +
  #geom_errorbar(size=0.8, width=0.2, aes(ymax=..y..,ymin=..y.., y = mean)) +
  theme(legend.position = "none")
#Diamonds are TX12, circles CA09
figure2b
ggsave(figure2b, filename = "figure2b.pdf")

#Add bioreplicate
model1_segment_dvg_all <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = gcay.screen.count.sgmt)
summary(model1_segment_dvg_all)
summary.aov(model1_segment_dvg_all)
#Multiple R-squared:  0.1409,	Adjusted R-squared:  0.1213 
#F-statistic: 7.176 on 16 and 700 DF,  p-value: 1.418e-15
#Significant, hits still there on Insu, MK2206, Ado and Uri. Alpe borderline

#Add strain
model2_segment_dvg_all <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = gcay.screen.count.sgmt)
summary(model2_segment_dvg_all)
summary.aov(model2_segment_dvg_all)
#Interesting that strain doesn't come up as significant and basically the same treatments come up as significant

(summary.aov(model2_segment_dvg_all)[[1]]$'Sum Sq' / sum(summary.aov(model2_segment_dvg_all)[[1]]$'Sum Sq'))*100
#13.11796719  0.97310998  0.05897013 85.84995270


#Comparisons 
model2_segment_dvg_all.glht.ctrl <- glht(model2_segment_dvg_all, linfct = mcp(Tx= "Dunnet"))
summary(model2_segment_dvg_all.glht.ctrl)
#Linear Hypotheses:
#Estimate Std. Error t value Pr(>|t|)  
#Alpe - DMSO == 0      32675      38086   0.858   0.8550  
#MK2206 - DMSO == 0   -18943      38086  -0.497   0.9802  
#4-OI - DMSO == 0     103797      38086   2.725   0.0769 .
#UK5099 - DMSO == 0    44272      38086   1.162   0.6732  
#Favp - DMSO == 0      60615      38086   1.592   0.4149 

#4-OI now drops out as borderline

#Now, simultaneous confidence intervals
confint(model1_TVG_overall_ca_dmso.glht.ctrl)
plot(confint(model1_TVG_overall_ca_dmso.glht.ctrl))

#Strain is NS, but how does a model with TX12 and CA09 separately look?
#CA09
model0_segment_dvg_ca <- lm(prop.dvg.sgmt ~ Tx, data = subset(gcay.screen.count.sgmt, host.ptg = "MDCK.CA09(H1N1)"))
summary(model0_segment_dvg_ca)
summary.aov(model0_segment_dvg_ca)
#Same hits as overall model

model0_segment_dvg_tx <- lm(prop.dvg.sgmt ~ Tx, data = subset(gcay.screen.count.sgmt, host.ptg = "MDCK.TX12(H3N2)"))
summary(model0_segment_dvg_tx)
summary.aov(model0_segment_dvg_tx)
#Same hits as overall model

#Same hits in both models tell me that there is a consistent effect of treatment across strains in DVGs

#Let's re-level segments:
levels(gcay.screen.count.sgmt$geno.sgmt)
#[1] "NS"  "PB2" "PB1" "PA"  "HA"  "NP"  "NA"  "M"
gcay.screen.count.sgmt$geno.sgmt <- fct_relevel(gcay.screen.count.sgmt$geno.sgmt, "PB2" , "PB1" , "PA" , "HA" , "NP" , "NA" , "M" , "NS")
levels(gcay.screen.count.sgmt$geno.sgmt)
gcay.screen.count.sgmt <- within(gcay.screen.count.sgmt, geno.sgmt <- relevel(geno.sgmt, ref = 8))

#Add segment
model3_segment_dvg_all <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg + geno.sgmt, data = gcay.screen.count.sgmt)
summary(model3_segment_dvg_all)
summary.aov(model3_segment_dvg_all)
#Segment is obviously significant. All the same results hold, but Alpe is now significant

#Add segment:Tx interaction
model4_segment_dvg_all <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg + geno.sgmt + geno.sgmt:Tx, data = gcay.screen.count.sgmt)
summary(model4_segment_dvg_all)
summary.aov(model4_segment_dvg_all)
#Now the overall effect of some of the individual treatments is gone, but some hits with Tx:geno.sgmt
#But need to control with each vehicle

#This is the "ideal" model

#Curious with what step will do
step(model4_segment_dvg_all)
#Only dropped host!


## Let's subset data. Multiple comparisons against a control should be done with each vehicle, but
  ## leaving the strain out for now 

#Let's split up data frame into host and vehicle
dvg_sgmt_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O")

dvg_sgmt_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO")

#Near-complete model for dH20
model4_segment_dvg_dH20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dh20)
summary(model4_segment_dvg_dH20)
summary.aov(model4_segment_dvg_dH20)
#Residual standard error: 0.1291 on 354 degrees of freedom
#Multiple R-squared:  0.6547,	Adjusted R-squared:  0.5825 
#F-statistic: 9.071 on 74 and 354 DF,  p-value: < 2.2e-16

#Significant, hits mostly on polymerase segments in insulin, ado and uri
#TxInsu:geno.sgmtPB2, TxAdo:geno.sgmtPB2, TxInsu:geno.sgmtPB1, TxAdo:geno.sgmtPB1, TxUri:geno.sgmtPB1, TxInsu:geno.sgmtPA, TxAdo:geno.sgmtPA. 
#Again, strain is not significant

#Near-complete model for dmso
model4_segment_dvg_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dmso)
summary(model4_segment_dvg_dmso)
summary.aov(model4_segment_dvg_dmso)
#Interesting, Tx:geno.sgmt is not significant and only get the TxMK2206:geno.sgmtPB1 as a hit and its a decrease

## Total Count overall 
model0_segment_tvg_all <- lm(count.tvg.sgmt ~ Tx, data = gcay.screen.count.sgmt)
summary(model0_segment_tvg_all)
summary.aov(model0_segment_tvg_all)

## Total Count by vehicle control for bioreplicate
model1_segment_tvg_all <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = gcay.screen.count.sgmt)
summary(model1_segment_tvg_all)
summary.aov(model1_segment_tvg_all)

## Total Count by vehicle control for bioreplicate and strain
model2_segment_tvg_all <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = gcay.screen.count.sgmt)
summary(model2_segment_tvg_all)
summary.aov(model2_segment_tvg_all)

(summary.aov(model2_segment_tvg_all)[[1]]$'Sum Sq' / sum(summary.aov(model2_segment_tvg_all)[[1]]$'Sum Sq'))*100
#8.0104326  1.1032324  0.7949817 90.0913534


#Strain is significant, so let's break down by strain

### Total Count by vehicle control for bioreplicate CA09
model2_segment_tvg_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = subset(gcay.screen.count.sgmt, host.ptg = "MDCK.CA09(H1N1)"))
summary(model2_segment_tvg_ca)
summary.aov(model2_segment_tvg_ca)

model2_segment_tvg_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = subset(gcay.screen.count.sgmt, host.ptg = "MDCK.TX12(H3N2)"))
summary(model2_segment_tvg_tx)
summary.aov(model2_segment_tvg_tx)

## Total Count by vehicle

# Total genome count Near-complete model for dH20
model4_segment_tvg_dH20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dh20)
summary(model4_segment_tvg_dH20)
summary.aov(model4_segment_tvg_dH20)
#Multiple R-squared:  0.5491,	Adjusted R-squared:  0.4548 
#F-statistic: 5.825 on 74 and 354 DF,  p-value: < 2.2e-16

#Significant and a bunch of interesting results
#Ado increases count
#First strain difference is that TX12 has lower total counts, so might make sense to split
#Now no Tx:geno.sgmt interactions, so this is an argument to split into strain as well which makes sense given different starting yields

#Near-complete model for dmso
model4_segment_tvg_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dmso)
summary(model4_segment_tvg_dmso)
summary.aov(model4_segment_tvg_dmso)
#Multiple R-squared:  0.5759,	Adjusted R-squared:  0.4864 
#F-statistic: 6.437 on 50 and 237 DF,  p-value: < 2.2e-16

#Significant and more interesting results 
#4-OI, UK5099, and Favp increase TVGs
#TX has lower total but borderline
#Hits on Tx4-OI:geno.sgmtPB2, Tx4-OI:geno.sgmtPB1, Tx4-OI:geno.sgmtPA all decreasing, but Tx:geno.sgmt not significant

## Total Count by vehicle and strain

#First CA09
#Let's split up data frame into host and vehicle
dvg_sgmt_dh20_ca <-  dvg_sgmt_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")

dvg_sgmt_dmso_ca <-  dvg_sgmt_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")

#The parameter estimates were wonky previously because using NS as the reference group for DVGs, which makes sense as NS has the lowest DVGs
# However, in TVGs NS is actually the highest value 
#So here we will re-level and set PB1 (lowest TVGs) to be the reference

#Let's re-level segments:
levels(gcay.screen.count.sgmt$geno.sgmt)
#[1] "NS"  "PB2" "PB1" "PA"  "HA"  "NP"  "NA"  "M"
gcay.screen.count.sgmt$geno.sgmt <- fct_relevel(gcay.screen.count.sgmt$geno.sgmt, "PB2" , "PB1" , "PA" , "HA" , "NP" , "NA" , "M" , "NS")
levels(gcay.screen.count.sgmt$geno.sgmt)
#[1] "PB2" "PB1" "PA"  "HA"  "NP"  "NA"  "M"   "NS" 
#Ahhh, better
#Relevel to PB1
gcay.screen.count.sgmt <- within(gcay.screen.count.sgmt, geno.sgmt <- relevel(geno.sgmt, ref = 2))

#Subset
dvg_sgmt_dh20_ca <-  gcay.screen.count.sgmt %>%
  filter(host.ptg == "MDCK.CA09(H1N1)") %>%
  filter(veh.Tx == "dH2O")


dvg_sgmt_dmso_ca <-  gcay.screen.count.sgmt %>%
  filter(host.ptg == "MDCK.CA09(H1N1)") %>%
  filter(veh.Tx == "DMSO")

#First dH20 
model4_segment_tvg_dH20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dh20_ca)
summary(model4_segment_tvg_dH20_ca)
summary.aov(model4_segment_tvg_dH20_ca)
#Multiple R-squared:  0.6665,	Adjusted R-squared:  0.4951 
#F-statistic: 3.888 on 73 and 142 DF,  p-value: 2.277e-1

#Insulin, Ado and Uri increase tvg, no Tx:geno.sgmt interaction

#Second DMSO
model4_segment_tvg_dmso_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dmso_ca)
summary(model4_segment_tvg_dmso_ca)
summary.aov(model4_segment_tvg_dmso_ca)
#Multiple R-squared:  0.7116,	Adjusted R-squared:  0.5613 
#F-statistic: 4.734 on 49 and 94 DF,  p-value: 5.367e-11

#Hits on 4-OI, UK5099, Favp increasing TVG
#These go away if I put PB1 as the reference

#Hits on Tx4-OI:geno.sgmtPB2, Tx4-OI:geno.sgmtPB1, Tx4-OI:geno.sgmtPA, Tx4-OI:geno.sgmtHA, Tx4-OI:geno.sgmtNP
#if I put PB1 as the reference now it's the mirror image: 
# Tx4-OI:geno.sgmtM      29168.7    10313.2   2.828  0.00572 **
#Tx4-OI:geno.sgmtNS     31431.0    10313.2   3.048  0.00299 **

#Second TX12
#Let's split up data frame into host and vehicle
dvg_sgmt_dh20_tx <-  dvg_sgmt_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)") %>%
  filter(veh.Tx == "dH2O")

dvg_sgmt_dmso_tx <-  dvg_sgmt_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)") %>%
  filter(veh.Tx == "dmso")

#First dH20 
model4_segment_tvg_dH20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dh20_tx)
summary(model4_segment_tvg_dH20_tx)
summary.aov(model4_segment_tvg_dH20_tx)
#Multiple R-squared:  0.7542,	Adjusted R-squared:  0.6251 
#F-statistic: 5.843 on 73 and 139 DF,  p-value: < 2.2e-16

#Model significant, but only hit TxTolyp:geno.sgmtM

#Second DMSO
model4_segment_tvg_dmso_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dmso_tx)
summary(model4_segment_tvg_dmso_tx)
summary.aov(model4_segment_tvg_dmso_tx)
#Multiple R-squared:  0.6444,	Adjusted R-squared:  0.459 
#F-statistic: 3.476 on 49 and 94 DF,  p-value: 1.083e-07

#Significant, but similarly few hits only  Tx4-OI:geno.sgmtM, TxUK5099:geno.sgmtM 

## Prop.DVG by vehicle and strain (by vehicle only above)

#Let's remember to re-level back to set NS as reference
#Let's re-level segments:
gcay.screen.count.sgmt$geno.sgmt <- fct_relevel(gcay.screen.count.sgmt$geno.sgmt, "PB2" , "PB1" , "PA" , "HA" , "NP" , "NA" , "M" , "NS")
levels(gcay.screen.count.sgmt$geno.sgmt)
gcay.screen.count.sgmt <- within(gcay.screen.count.sgmt, geno.sgmt <- relevel(geno.sgmt, ref = 8))

#First CA09
#First dH20 CA09 PropDVG
model4_segment_dvg_dH20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dh20_ca)
summary(model4_segment_dvg_dH20_ca)
summary.aov(model4_segment_dvg_dH20_ca)
#Multiple R-squared:  0.6583,	Adjusted R-squared:  0.4826 
#F-statistic: 3.747 on 73 and 142 DF,  p-value: 8.368e-12

#Hits increases in TxAdo:geno.sgmtPB2, TxInsu:geno.sgmtPB1, TxAdo:geno.sgmtPB1, TxSynec:geno.sgmtPB1, TxAdo:geno.sgmtPA, TxAdo:geno.sgmtHA

#Second DMSO CA09 PropDVG
model4_segment_dvg_dmso_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dmso_ca)
summary(model4_segment_dvg_dmso_ca)
summary.aov(model4_segment_dvg_dmso_ca)
#NS

#Second TX12
#First TX12 dH20 Prop DVG
model4_segment_dvg_dH20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dh20_tx)
summary(model4_segment_dvg_dH20_tx)
summary.aov(model4_segment_dvg_dH20_tx)
#Multiple R-squared:  0.7766,	Adjusted R-squared:  0.6593 
#F-statistic:  6.62 on 73 and 139 DF,  p-value: < 2.2e-16

#Model significant, hits for increases on TxInsu:geno.sgmtPB2, TxAdo:geno.sgmtPB2, TxInsu:geno.sgmtPB1, TxAdo:geno.sgmtPB1, TxInsu:geno.sgmtPA, TxAdo:geno.sgmtPA 
#These increases are also just very big

#Second TX12 DMSO Prop DVG 
model4_segment_dvg_dmso_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + geno.sgmt + geno.sgmt:Tx, data = dvg_sgmt_dmso_tx)
summary(model4_segment_dvg_dmso_tx)
summary.aov(model4_segment_dvg_dmso_tx)
#Multiple R-squared:  0.4538,	Adjusted R-squared:  0.1691 
#F-statistic: 1.594 on 49 and 94 DF,  p-value: 0.02691

#Significant only hit increase in TxMK2206:geno.sgmtPB1


##### Overall Summary so far
#Looks like TVG are a bust using the overall data set for both strains
#prop DVG provided some hits for Texas, namely Insulin, MK2206, Ado, Uri but Insulin and MK2206 drop off as borderline when controlling for vehicle

#Seems to me that it's better to use the per-segment data set in any case 
#Per-segment PropDVG
#Using per segment data set we have quite a few significant results and things to note
#For propDVG Insu, MK2206, Ado and Uri (Alpe inconsistent-ish) came out without controlling for much at all and across strains
  #Super noteworthy that strain at this point is not significant at all
#Then we drill down for specific comparisons to test the effect of each one, so splitting by vehicle
# notably not splitting by strain here
#For H20: Significant, hits on polymerase segments in insulin, ado and uri
#TxInsu:geno.sgmtPB2, TxAdo:geno.sgmtPB2, TxInsu:geno.sgmtPB1, TxAdo:geno.sgmtPB1, TxUri:geno.sgmtPB1, TxInsu:geno.sgmtPA, TxAdo:geno.sgmtPA. 
#For DMSO: TxMK2206:geno.sgmtPB1
#Splitting by strain as a double check, but NS in overall model

#Per-segment Total VG

#Here splitting by vehicle and strain because TVG levels are so different and vehicle needed for comparisons

#CA09 dH20: Insulin, Ado and Uri increase tvg, no Tx:geno.sgmt interaction
#CA09 DMSO:  4-OI, UK5099, Favp increasing TVG, also #Hits on Tx4-OI:geno.sgmtPB2, Tx4-OI:geno.sgmtPB1, Tx4-OI:geno.sgmtPA, Tx4-OI:geno.sgmtHA, Tx4-OI:geno.sgmtNP

#TX12 dH20: Model significant, but only hit TxTolyp:geno.sgmtM
#TX12 DMSO: Significant, but similarly few hits only  Tx4-OI:geno.sgmtM, TxUK5099:geno.sgmtM 

#But then it's difficult to compare TVG and prop:DVG

#Here's propDVG Results

#CA09 dH20: Significant, Hits increases in TxAdo:geno.sgmtPB2, TxInsu:geno.sgmtPB1, TxAdo:geno.sgmtPB1, TxSynec:geno.sgmtPB1, TxAdo:geno.sgmtPA, TxAdo:geno.sgmtHA
#CA09 DMSO: Model NS and no hits

#TX12 DH20 Model significant, hits for increases on TxInsu:geno.sgmtPB2, TxAdo:geno.sgmtPB2, TxInsu:geno.sgmtPB1, TxAdo:geno.sgmtPB1, TxInsu:geno.sgmtPA, TxAdo:geno.sgmtPA These increases are also just very big
#TX12 DMSO Significant only hit increase in TxMK2206:geno.sgmtPB1

#The only significant non-polymerase segment result was prop DVG TxAdo:geno.sgmtHA in CA09 and
#TVG Tx4-OI:geno.sgmtHA, Tx4-OI:geno.sgmtNP for CA09 and Tx4-OI:geno.sgmtM, TxUK5099:geno.sgmtM with TX12

#Might want to do some testing of Alpelisib on its own to verify results from our paper

###### / Stat Testing ######


##### Stat Plots ####

#Plots that reflect the analysis

#Another view not accounting for strain Proportion of DVGs regardless of strain / vehicle showing the players
#Players by eye are Insulin, Alpe MK2206, 4-OI, UK5099, Ado, Uri,
ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate, shape=host.ptg)) +
  geom_point() + theme(legend.position = "none")

#Now 
ggplot(gcay.screen.count.sgmt, aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate, shape=host.ptg)) +
  geom_point() + theme(legend.position = "none")

#Prop DVG Polymerase segments only by vehicle and strain
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#Prop DVG Other segments by vehicle and strain
non_polymerase_segments <- c("HA", "NP", "NA", "M", "NS")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  non_polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#TVG Polymerase segments only by vehicle and strain
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#TVG Other segments by vehicle and strain
non_polymerase_segments <- c("HA", "NP", "NA", "M", "NS")
gcay.screen.count.sgmt %>%
  filter(geno.sgmt %in%  non_polymerase_segments) %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")


#Adenosine and insulin
#Vehicle combinations propDelVG
ado_insu <- c("dH2O", "Ado", "Insu")
polymerase_segments <- c("PB2", "PB1", "PA")
figure3 <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  ado_insu) %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  mutate(Tx = fct_relevel(Tx, c("dH2O", "Ado", "Insu")))  %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  ylim(0, 1.0) +
  xlab("Viral Genome Segment") +
  ylab("Proportion of Defective Viral Genomes (DelVGs)") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  theme(panel.margin = unit(3, "lines")) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + theme(legend.position = "none") +
  scale_fill_colorblind() +
  facet_grid(host.ptg~Tx, scales = "free")

figure3
ggsave(figure3, filename = "figure3.pdf")

#Gel analysis
#Followup gel experiment looking at wild-type TX12 infections under increasing doses of adenosine
gelgenie_qupath_band_data <- read.csv("data/gelgenie_qupath_band_data.csv")
colnames(gelgenie_qupath_band_data)[1] <- "Lane"
#Adding Sample names manually
Sample.Names <- c("TX12 WT No Ado", "TX12 WT 1µM Ado", "TX12 WT 5µM Ado", "TX12 WT 10µM Ado", "TX12 WT 20µM Ado")

#Binding sample names
cbind(Sample.Names, gelgenie_qupath_band_data)

#Local Corrected Volume
figure_Sgela <- ggplot(gelgenie_qupath_band_data, aes(x = fct_reorder(Sample.Names, Lane), y=Local.Corrected.Volume)) +
  ggtitle("Local Corrected Volume") +
  ylab("Intensity (A.U)") +
  xlab("Concentration of Adenosine (Ado) in WT TX12 Infections") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  theme(panel.margin = unit(3, "lines")) +
  geom_col()

#Rolling Ball Corrected Volume
figure_Sgelb <- ggplot(gelgenie_qupath_band_data, aes(x = fct_reorder(Sample.Names, Lane), y=Rolling.Ball.Corrected.Volume)) +
  ggtitle("Rolling Ball Corrected Volume") +
  ylab("Intensity (A.U)") +
  xlab("Concentration of Adenosine (Ado) in WT TX12 Infections") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  theme(panel.margin = unit(3, "lines")) +
  geom_col()

# 4-OI, UK5099, and Favi - TVG 
tvg_increasers <- c("DMSO", "4-OI", "UK5099", "Favp")
figure4 <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  tvg_increasers) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  xlab("Viral Genome Segment") +
  ylab("Total Viral Genomes") +
  ylim(0, 80000) +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  theme(panel.margin = unit(3, "lines")) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + scale_fill_colorblind() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")
figure4
ggsave(figure4, filename = "figure4.pdf")

# 4-OI, UK5099, and Favi - TVG 
cyano <- c("dH2O", "Lepto", "Nosto", "Oscil", "Synec", "Tolyp")
antigenic_segments <- c("HA", "NA")
figure5 <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  cyano) %>%
  filter(geno.sgmt %in% antigenic_segments) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  ylim(0, 15000) +
  xlab("Viral Genome Segment") +
  ylab("Total Viral Genomes") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  theme(panel.margin = unit(3, "lines")) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") +
  scale_fill_colorblind() +
  theme(legend.position = "none") + 
  facet_grid(host.ptg~Tx, scales = "free")
figure5
ggsave(figure5, filename = "figure5.pdf")

##### / Stat Plots ####

##### Control Analyses and sequencing stats ####
demultiplexing_trimming_stats <- read.csv("~/Dropbox/mixtup/Documentos/ucdavis/papers/dvg_drug_screen/influenza_dvg_drug_screen/output/demultiplexing_trimming_stats.csv", header=FALSE)

#Quick and dirty removing negative controls "by hand"
demultiplexing_trimming_stats[c(-44, -56, -48, -40, -52, -60, -97), ]$V6
mean(demultiplexing_trimming_stats[c(-44, -56, -48, -40, -52, -60, -97), ]$V6)
#[1] 67496.81

sd(demultiplexing_trimming_stats[-c(44, 56, 48, 40, 52, 60, 97), ]$V6)
#56183.6

## Control Plot
vehicles <- c("dH2O", "DMSO")

#All four strain/vehicle combinations propDelVG
gcay.screen.count.sgmt %>%
  filter(Tx %in%  vehicles) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_wrap(host.ptg~Tx, scales = "free")

#All four strain/vehicle combinations TVG
gcay.screen.count.sgmt %>%
  filter(Tx %in%  vehicles) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_wrap(host.ptg~Tx, scales = "free")

#Colorblind friendly palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
#and http://jfly.iam.u-tokyo.ac.jp/color/:

#NOV 14 START HERE
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)

#Pick CA09 as representative
#Vehicle combinations propDelVG
figure1a <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  vehicles) %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")  %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  ylim(0, 1.0) +
  xlab("Viral Genome Segment") +
  ylab("Proportion of Defective Viral Genomes (DelVGs)") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + scale_fill_colorblind() + theme(legend.position = "none") + facet_grid(.~Tx, scales = "free")
figure1a
ggsave(figure1a, filename = "figure1a.pdf")

#Quick for Ile
mk2206 <- c("DMSO", "MK2206")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  mk2206) %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")  %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate)) +
  #ylim(0, 1.0) +
  xlab("Viral Genome Segment") +
  #ylab("Proportion of Defective Viral Genomes (DelVGs)") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + theme(legend.position = "none") + facet_grid(geno.sgmt~., scales = "free")



#Pick CA09 as representative
#vehicle combinations TVG
figure1b <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  vehicles) %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")  %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  ylim(0, 80000) +
  xlab("Viral Genome Segment") +
  ylab("Total Viral Genomes") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + scale_fill_colorblind() + theme(legend.position = "none") + facet_grid(.~Tx, scales = "free")
figure1b
ggsave(figure1b, filename = "figure1b.pdf")

#Creating a supplement of Figure 1 which includes the TX12 data
#Vehicle combinations propDelVG
figureS6a <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  vehicles) %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")  %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=prop.dvg.sgmt, color=bioreplicate)) +
  ylim(0, 1.0) +
  xlab("Viral Genome Segment") +
  ylab("Proportion of Defective Viral Genomes (DelVGs)") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + scale_fill_colorblind() + theme(legend.position = "none") + facet_grid(.~Tx, scales = "free")
figureS6a
ggsave(figureS6a, filename = "figureS6a.pdf")
#Second panel for TX12 Supplement
#vehicle combinations TVG
figureS6b <- gcay.screen.count.sgmt %>%
  filter(Tx %in%  vehicles) %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")  %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  ylim(0, 80000) +
  xlab("Viral Genome Segment") +
  ylab("Total Viral Genomes") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + scale_fill_colorblind() + theme(legend.position = "none") + facet_grid(.~Tx, scales = "free")
figureS6b
ggsave(figureS6b, filename = "figureS6b.pdf")

#MK2206 Plot Polymerase
mk2206 <- c("DMSO", "MK2206")
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  mk2206) %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(geno.sgmt~host.ptg, scales = "free")

#MK2206 Plot Non-Polymerase
mk2206 <- c("DMSO", "MK2206")
non_polymerase_segments <- c("HA", "NP", "NA", "M", "NS")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  mk2206) %>%
  filter(geno.sgmt %in%  non_polymerase_segments) %>%
  ggplot(aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(geno.sgmt~host.ptg, scales = "free")


#Adenosine Plot Polymerase
ado <- c("DMSO", "Ado")
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  ado) %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(geno.sgmt~host.ptg, scales = "free")

#Adenosine Plot Non-Polymerase
ado <- c("DMSO", "Ado")
non_polymerase_segments <- c("HA", "NP", "NA", "M", "NS")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  ado) %>%
  filter(geno.sgmt %in%  non_polymerase_segments) %>%
  ggplot(aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(geno.sgmt~host.ptg, scales = "free")


#4-OI Plot Polymerase TVG
four_oi <- c("DMSO", "4-OI")
polymerase_segments <- c("PB2", "PB1", "PA")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  four_oi) %>%
  filter(geno.sgmt %in%  polymerase_segments) %>%
  ggplot(aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(geno.sgmt~host.ptg, scales = "free")


#4-OI Plot TVG Per Segment
four_oi <- c("DMSO", "4-OI")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  four_oi) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#4-OI Plot TVG All Segments Together
four_oi <- c("DMSO", "4-OI")
gcay.screen.count.sgmt %>%
  filter(Tx %in%  four_oi) %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=Tx, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~., scales = "free")

#TVG All Segments Together All Treatments
gcay.screen.count.sgmt %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")

#####  / Control Analyses and sequencing stats ####



#Shit should I have been making models by *segment* instead?

### First let's do water ####
## PB2
dvg_pb2_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "PB2")

model0_dvg_pb2_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_pb2_dh20)
summary(model0_dvg_pb2_dh20)
summary.aov(model0_dvg_pb2_dh20)
#TxInsu                   0.242260   0.080745   3.000  0.00452 ** 
#  TxAdo                    0.476894   0.080745   5.906 5.42e-07 ***
#  TxUri                    0.151488   0.080745   1.876  0.06760 .
#host.ptgMDCK.TX12(H3N2)  0.081063   0.038063   2.130  0.03910 *

# CA09
dvg_pb2_dh20_ca <- dvg_pb2_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")

model0_dvg_pb2_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_pb2_dh20_ca)
summary(model0_dvg_pb2_dh20_ca)
summary.aov(model0_dvg_pb2_dh20_ca)
#TxAdo          0.350356   0.107370   3.263  0.00488 **
#But overall model NS
#F-statistic: 2.037 on 10 and 16 DF,  p-value: 0.09861

#TX12
dvg_pb2_dh20_tx <- dvg_pb2_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")

model0_dvg_pb2_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_pb2_dh20_tx)
summary(model0_dvg_pb2_dh20_tx)
summary.aov(model0_dvg_pb2_dh20_tx)
#TxInsu         0.365633   0.122571   2.983 0.008785 ** 
#  TxAdo          0.603432   0.122571   4.923 0.000153 ***

## PB1
dvg_pb1_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "PB1")

model0_dvg_pb1_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_pb1_dh20)
summary(model0_dvg_pb1_dh20)
summary.aov(model0_dvg_pb1_dh20)
#TxInsu                   0.424720   0.137739   3.084  0.00361 ** 
#  TxAdo                    0.692711   0.137739   5.029 9.68e-06 ***
#  TxUri                    0.233630   0.137739   1.696  0.09725 . 

# CA09
dvg_pb1_dh20_ca <- dvg_pb1_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")

model0_dvg_pb1_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_pb1_dh20_ca)
summary(model0_dvg_pb1_dh20_ca)
summary.aov(model0_dvg_pb1_dh20_ca)
#TxAdo          0.575611   0.204348   2.817   0.0124 *


#TX12
dvg_pb1_dh20_tx <- dvg_pb1_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")

model0_dvg_pb1_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_pb1_dh20_tx)
summary(model0_dvg_pb1_dh20_tx)
summary.aov(model0_dvg_pb1_dh20_tx)
#TxInsu         0.52964    0.20321   2.606  0.01909 * 
#TxAdo          0.80981    0.20321   3.985  0.00106 *


## PA
dvg_pa_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "PA")

model0_dvg_pa_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_pa_dh20)
summary(model0_dvg_pa_dh20)
summary.aov(model0_dvg_pa_dh20)
#TxInsu                   0.255170   0.109453   2.331   0.0246 *  
#  TxAdo                    0.650365   0.109453   5.942 4.81e-07 ***
#  TxUri                    0.185146   0.109453   1.692   0.0981 .

# CA09
dvg_pa_dh20_ca <- dvg_pa_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_pa_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_pa_dh20_ca)
summary(model0_dvg_pa_dh20_ca)
summary.aov(model0_dvg_pa_dh20_ca)
#TxAdo          0.514519   0.140622   3.659  0.00212 **


#TX12
dvg_pa_dh20_tx <- dvg_pa_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")

model0_dvg_pa_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_pa_dh20_tx)
summary(model0_dvg_pa_dh20_tx)
summary.aov(model0_dvg_pa_dh20_tx)
#TxInsu         0.387404   0.178854   2.166 0.045756 *  
#TxAdo          0.786210   0.178854   4.396 0.000451 ***

## HA
dvg_ha_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "HA")

model0_dvg_ha_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_ha_dh20)
summary(model0_dvg_ha_dh20)
summary.aov(model0_dvg_ha_dh20)
#TxAdo                    0.184594   0.066058   2.794 0.007871 **
#host.ptgMDCK.TX12(H3N2) -0.116493   0.031509  -3.697 0.000639 ***

# CA09
dvg_ha_dh20_ca <- dvg_ha_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_ha_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_ha_dh20_ca)
summary(model0_dvg_ha_dh20_ca)
summary.aov(model0_dvg_ha_dh20_ca)
#TxAdo          0.3372571  0.1237154   2.726    0.015 *

#TX12
dvg_ha_dh20_tx <- dvg_ha_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")

model0_dvg_ha_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_ha_dh20_tx)
summary(model0_dvg_ha_dh20_tx)
summary.aov(model0_dvg_ha_dh20_tx)
#TxUri          0.089077   0.033944   2.624   0.0192 *
#But model NS

## NP
dvg_np_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "NP")

model0_dvg_np_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_np_dh20)
summary(model0_dvg_np_dh20)
summary.aov(model0_dvg_np_dh20)
#TxAdo                    0.035771   0.014293   2.503   0.0163 * 
#host.ptgMDCK.TX12(H3N2) -0.055049   0.006738  -8.170 3.24e-10 ***

#CA09
dvg_np_dh20_ca <- dvg_np_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_np_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_np_dh20_ca)
summary(model0_dvg_np_dh20_ca)
summary.aov(model0_dvg_np_dh20_ca)
#TxAdo          0.061575   0.028392   2.169   0.0455 *
#But model overall NS

#TX12
dvg_np_dh20_tx <- dvg_np_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_np_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_np_dh20_tx)
summary(model0_dvg_np_dh20_tx)
summary.aov(model0_dvg_np_dh20_tx)
#TxInsu         0.0077981  0.0036557   2.133 0.048746 *  
#  TxAdo          0.0099666  0.0036557   2.726 0.014943 *

## NA
dvg_na_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "NA")

model0_dvg_na_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_na_dh20)
summary(model0_dvg_na_dh20)
summary.aov(model0_dvg_na_dh20)
#TxLepto                 -9.685e-03  3.378e-03  -2.867  0.00658 **
#TxTolyp                 -5.896e-03  3.378e-03  -1.745  0.08862 . 
#host.ptgMDCK.TX12(H3N2)  9.472e-03  1.640e-03   5.776 9.77e-07 ***

#CA09
dvg_na_dh20_ca <- dvg_na_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_na_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_na_dh20_ca)
summary(model0_dvg_na_dh20_ca)
summary.aov(model0_dvg_na_dh20_ca)
#Nada

#TX12
dvg_na_dh20_tx <- dvg_na_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_na_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_na_dh20_tx)
summary(model0_dvg_na_dh20_tx)
summary.aov(model0_dvg_na_dh20_tx)
#TxLepto       -0.0194644  0.0047130  -4.130 0.001021 **
#TxTolyp       -0.0115042  0.0047130  -2.441 0.028536 * 

## M
dvg_m_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "M")

model0_dvg_m_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_m_dh20)
summary(model0_dvg_m_dh20)
summary.aov(model0_dvg_m_dh20)
#TxInsu                   1.401e-03  8.128e-04   1.724 0.092118 .
#TxUri                    1.569e-03  8.128e-04   1.931 0.060283 .  
#TxLepto                  1.535e-03  8.128e-04   1.888 0.065907 .
#host.ptgMDCK.TX12(H3N2)  1.242e-02  3.831e-04  32.415  < 2e-16 ***

#CA09
dvg_m_dh20_ca <- dvg_m_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_m_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_m_dh20_ca)
summary(model0_dvg_m_dh20_ca)
summary.aov(model0_dvg_m_dh20_ca)
#TxNosto       -7.990e-04  3.960e-04  -2.018   0.0607 .

#TX12
dvg_m_dh20_tx <- dvg_m_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_m_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_m_dh20_tx)
summary(model0_dvg_m_dh20_tx)
summary.aov(model0_dvg_m_dh20_tx)
#TxAdo         0.0027815  0.0015025   1.851   0.0827 .  
#TxUri         0.0029539  0.0015025   1.966   0.0669 .  
#TxLepto       0.0031102  0.0015025   2.070   0.0550 .  
#TxNosto       0.0029586  0.0015025   1.969   0.0665 .
#Overall model also NS

## NS
dvg_ns_dh20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "NS")

model0_dvg_ns_dh20 <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_ns_dh20)
summary(model0_dvg_ns_dh20)
summary.aov(model0_dvg_ns_dh20)
#TxLepto                  0.0019610  0.0011388   1.722   0.0924 .  
#TxNosto                  0.0027982  0.0011388   2.457   0.0182 *  
#  TxOscil                  0.0024397  0.0011388   2.142   0.0380 *  
#  TxSynec                  0.0028355  0.0011388   2.490   0.0168 *
#host.ptgMDCK.TX12(H3N2)  0.0141545  0.0005368  26.368   <2e-16 ***

#CA09
dvg_ns_dh20_ca <- dvg_ns_dh20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_ns_dh20_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_ns_dh20_ca)
summary(model0_dvg_ns_dh20_ca)
summary.aov(model0_dvg_ns_dh20_ca)
#Nada

#TX12
dvg_ns_dh20_tx <- dvg_ns_dh20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_ns_dh20_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_ns_dh20_tx)
summary(model0_dvg_ns_dh20_tx)
summary.aov(model0_dvg_ns_dh20_tx)
#TxLepto        0.0037638  0.0020031   1.879   0.0786 .  
#TxNosto        0.0052008  0.0020031   2.596   0.0195 *  
# TxOscil        0.0043773  0.0020031   2.185   0.0441 *  
# TxSynec        0.0056071  0.0020031   2.799   0.0129 *
#Overall model just borderline though

#Let's try a plot here:
dvg_ns_dh20 %>%
  ggplot(aes(x=Tx, y=prop.dvg.sgmt, color=bioreplicate)) +
  geom_point() + theme(legend.position = "none") + facet_grid(host.ptg~., scales = "free")
  

### Now let's do DMSO ####
## PB2
dvg_pb2_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "PB2")

model0_dvg_pb2_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_pb2_dmso)
summary(model0_dvg_pb2_dmso)
summary.aov(model0_dvg_pb2_dmso)
#Nada

## PB1
dvg_pb1_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "PB1")

model0_dvg_pb1_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_pb1_dmso)
summary(model0_dvg_pb1_dmso)
summary.aov(model0_dvg_pb1_dmso)
#Nada

## PA
dvg_pa_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "PA")

model0_dvg_pa_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_pa_dmso)
summary(model0_dvg_pa_dmso)
summary.aov(model0_dvg_pa_dmso)
#Nada

## HA
dvg_ha_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "HA")

model0_dvg_ha_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_ha_dmso)
summary(model0_dvg_ha_dmso)
summary.aov(model0_dvg_ha_dmso)
#nada

## NP
dvg_np_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "NP")

model0_dvg_np_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_np_dmso)
summary(model0_dvg_np_dmso)
summary.aov(model0_dvg_np_dmso)
#host.ptgMDCK.TX12(H3N2) -0.048634   0.006156  -7.900 1.71e-08 ***
#No segment specific effects

#CA09
dvg_np_dmso_ca <- dvg_np_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_np_dmso_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_np_dmso_ca)
summary(model0_dvg_np_dmso_ca)
summary.aov(model0_dvg_np_dmso_ca)
#Nada

#TX12
dvg_np_dmso_tx <- dvg_np_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_np_dmso_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_np_dmso_tx)
summary(model0_dvg_np_dmso_tx)
summary.aov(model0_dvg_np_dmso_tx)
#Nada

## NA
dvg_na_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "NA")

model0_dvg_na_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_na_dmso)
summary(model0_dvg_na_dmso)
summary.aov(model0_dvg_na_dmso)
#TxMK2206                 0.0074195  0.0034901   2.126   0.0428 * 
#host.ptgMDCK.TX12(H3N2)  0.0105239  0.0020150   5.223 1.67e-05 ***

#CA09
dvg_na_dmso_ca <- dvg_na_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_na_dmso_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_na_dmso_ca)
summary(model0_dvg_na_dmso_ca)
summary.aov(model0_dvg_na_dmso_ca)
#Nada

#TX12
dvg_na_dmso_tx <- dvg_na_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_na_dmso_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_na_dmso_tx)
summary(model0_dvg_na_dmso_tx)
summary.aov(model0_dvg_na_dmso_tx)
#TxMK2206       0.0142625  0.0047605   2.996   0.0134 * 
#TxFavp        -0.0121656  0.0047605  -2.556   0.0286 * 


## M
dvg_m_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "M")

model0_dvg_m_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_m_dmso)
summary(model0_dvg_m_dmso)
summary.aov(model0_dvg_m_dmso)
#host.ptgMDCK.TX12(H3N2)  0.0120192  0.0008294  14.491 2.96e-14 ***

#CA09
dvg_m_dmso_ca <- dvg_m_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_m_dmso_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_m_dmso_ca)
summary(model0_dvg_m_dmso_ca)
summary.aov(model0_dvg_m_dmso_ca)
#Nada / some borderline

#TX12
dvg_m_dmso_tx <- dvg_m_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_m_dmso_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_m_dmso_tx)
summary(model0_dvg_m_dmso_tx)
summary.aov(model0_dvg_m_dmso_tx)
#Nada

## NS
dvg_ns_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "NS")

model0_dvg_ns_dmso <- lm(prop.dvg.sgmt ~ Tx + bioreplicate + host.ptg, data = dvg_ns_dmso)
summary(model0_dvg_ns_dmso)
summary.aov(model0_dvg_ns_dmso)
#TxFavp                  -0.0029547  0.0016693  -1.770   0.0880 . 
#host.ptgMDCK.TX12(H3N2)  0.0139332  0.0009638  14.457 3.13e-14 ***

#CA09
dvg_ns_dmso_ca <- dvg_ns_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_dvg_ns_dmso_ca <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_ns_dmso_ca)
summary(model0_dvg_ns_dmso_ca)
summary.aov(model0_dvg_ns_dmso_ca)
#Nada

#TX12
dvg_ns_dmso_tx <- dvg_ns_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_dvg_ns_dmso_tx <- lm(prop.dvg.sgmt ~ Tx + bioreplicate, data = dvg_ns_dmso_tx)
summary(model0_dvg_ns_dmso_tx)
summary.aov(model0_dvg_ns_dmso_tx)
#TxFavp        -0.006108   0.002764  -2.210   0.0516 .

### TVG DMSO ####
## PB2
tvg_pb2_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "PB2")

model0_tvg_pb2_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_pb2_dmso)
summary(model0_tvg_pb2_dmso)
summary.aov(model0_tvg_pb2_dmso)
#Tx4-OI                    2287.3      905.6   2.526   0.0177 *
#TxFavp                    1576.7      905.6   1.741   0.0931 .


## PB1
tvg_pb1_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "PB1")

model0_tvg_pb1_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_pb1_dmso)
summary(model0_tvg_pb1_dmso)
summary.aov(model0_tvg_pb1_dmso)
#Tx4-OI                   1177.67     462.20   2.548   0.0168 *
#TxFavp                    934.67     462.20   2.022   0.0532 .


## PA
tvg_pa_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "PA")

model0_tvg_pa_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_pa_dmso)
summary(model0_tvg_pa_dmso)
summary.aov(model0_tvg_pa_dmso)
#Tx4-OI                    3886.3     1285.5   3.023  0.00543 **
#TxUK5099                  2610.3     1285.5   2.031  0.05225 . 
#TxFavp                    2767.7     1285.5   2.153  0.04041 *

## HA
tvg_ha_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "HA")

model0_tvg_ha_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_ha_dmso)
summary(model0_tvg_ha_dmso)
summary.aov(model0_tvg_ha_dmso)
#Tx4-OI                   8590.00    2447.24   3.510  0.00159 **
#TxUK5099                 5051.67    2447.24   2.064  0.04873 *

## NP
tvg_np_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "NP")

model0_tvg_np_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_np_dmso)
summary(model0_tvg_np_dmso)
summary.aov(model0_tvg_np_dmso)
#Tx4-OI                   9482.83    2792.07   3.396  0.00213 **
#TxUK5099                 6189.17    2792.07   2.217  0.03526 * 
#TxFavp                   5809.50    2792.07   2.081  0.04708 *

#These should be interpreted as the average of the increases / decreases

#CA09
tvg_np_dmso_ca <- tvg_np_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_np_dmso_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_np_dmso_ca)
summary(model0_tvg_np_dmso_ca)
summary.aov(model0_tvg_np_dmso_ca)
#

#TX12
tvg_np_dmso_tx <- tvg_np_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_np_dmso_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_np_dmso_tx)
summary(model0_tvg_np_dmso_tx)
summary.aov(model0_tvg_np_dmso_tx)

## NA
tvg_na_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "NA")

model0_tvg_na_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_na_dmso)
summary(model0_tvg_na_dmso)
summary.aov(model0_tvg_na_dmso)
#Tx4-OI                   11798.3     3463.2   3.407  0.00207 **
#TxUK5099                  7577.3     3463.2   2.188  0.03750 *

## M
tvg_m_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "M")

model0_tvg_m_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_m_dmso)
summary(model0_tvg_m_dmso)
summary.aov(model0_tvg_m_dmso)
#Tx4-OI                     30476       8857   3.441   0.0019 **
#TxUK5099                   20742       8857   2.342   0.0268 * 
#TxFavp                     19298       8857   2.179   0.0382 * 


## NS
tvg_ns_dmso <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "DMSO") %>%
  filter(geno.sgmt == "NS")

model0_tvg_ns_dmso <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_ns_dmso)
summary(model0_tvg_ns_dmso)
summary.aov(model0_tvg_ns_dmso)
#Tx4-OI                     20732       7340   2.824 0.008796 **
#TxFavp                     12762       7340   1.739 0.093515 .
#host.ptgMDCK.TX12(H3N2)   -15911       4238  -3.754 0.000845 ***

#CA09
tvg_ns_dmso_ca <- tvg_ns_dmso %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_ns_dmso_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_ns_dmso_ca)
summary(model0_tvg_ns_dmso_ca)
summary.aov(model0_tvg_ns_dmso_ca)
#Tx4-OI           33209      12276   2.705   0.0221 *

#TX12
tvg_ns_dmso_tx <- tvg_ns_dmso %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_ns_dmso_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_ns_dmso_tx)
summary(model0_tvg_ns_dmso_tx)
summary.aov(model0_tvg_ns_dmso_tx)
#Tx4-OI          8255.7     4169.2   1.980   0.0759 .
#TxUK5099        7932.7     4169.2   1.903   0.0862 .
#Overall model NS


### TVG H20 ####
## PB2
tvg_pb2_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "PB2")

model0_tvg_pb2_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_pb2_h20)
summary(model0_tvg_pb2_h20)
summary.aov(model0_tvg_pb2_h20)
#Nada

## PB1
tvg_pb1_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "PB1")

model0_tvg_pb1_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_pb1_h20)
summary(model0_tvg_pb1_h20)
summary.aov(model0_tvg_pb1_h20)
#TxTolyp                   701.17     352.36   1.990   0.0531 .
#host.ptgMDCK.TX12(H3N2)   424.44     166.10   2.555   0.0143 *

#CA09
tvg_pb1_h20_ca <- tvg_pb1_h20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_pb1_h20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_pb1_h20_ca)
summary(model0_tvg_pb1_h20_ca)
summary.aov(model0_tvg_pb1_h20_ca)
#Nada

#TX12
tvg_pb1_h20_tx <- tvg_pb1_h20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_pb1_h20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_pb1_h20_tx)
summary(model0_tvg_pb1_h20_tx)
summary.aov(model0_tvg_pb1_h20_tx)
#TxTolyp        1410.00     532.93   2.646   0.0176 *

## PA
tvg_pa_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "PA")

model0_tvg_pa_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_pa_h20)
summary(model0_tvg_pa_h20)
summary.aov(model0_tvg_pa_h20)
#TxTolyp                  1659.500    853.354   1.945  0.05853 .
#host.ptgMDCK.TX12(H3N2)   695.704    402.275   1.729  0.09108 .

#CA09
tvg_pa_h20_ca <- tvg_pa_h20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_pa_h20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_pa_h20_ca)
summary(model0_tvg_pa_h20_ca)
summary.aov(model0_tvg_pa_h20_ca)
#Nada

#TX12
tvg_pa_h20_tx <- tvg_pa_h20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_pa_h20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_pa_h20_tx)
summary(model0_tvg_pa_h20_tx)
summary.aov(model0_tvg_pa_h20_tx)
#TxTolyp         3444.0     1018.8   3.381  0.00381 **
#But the overall model wasn't significant to begin with

## HA
tvg_ha_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "HA")

model0_tvg_ha_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_ha_h20)
summary(model0_tvg_ha_h20)
summary.aov(model0_tvg_ha_h20)
#TxLepto                 -3029.50    1301.59  -2.328    0.025 *
#host.ptgMDCK.TX12(H3N2) -1290.12     620.84  -2.078    0.044 *

#CA09
tvg_ha_h20_ca <- tvg_ha_h20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_ha_h20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_ha_h20_ca)
summary(model0_tvg_ha_h20_ca)
summary.aov(model0_tvg_ha_h20_ca)
#Nada

#TX12
tvg_ha_h20_tx <- tvg_ha_h20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_ha_h20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_ha_h20_tx)
summary(model0_tvg_ha_h20_tx)
summary.aov(model0_tvg_ha_h20_tx)
#TxAdo         -2101.67    1052.26  -1.997 0.064266 .  
#TxUri         -4107.96    1191.07  -3.449 0.003580 ** 
#TxLepto       -4135.67    1052.26  -3.930 0.001336 **
#TxTolyp       -3178.00    1052.26  -3.020 0.008611 **

## NP
tvg_np_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "NP")

model0_tvg_np_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_np_h20)
summary(model0_tvg_np_h20)
summary.aov(model0_tvg_np_h20)
#TxTolyp                   3856.2     2048.4   1.883   0.0667 .
#host.ptgMDCK.TX12(H3N2)   1856.4      965.6   1.923   0.0613 .

#CA09
tvg_np_h20_ca <- tvg_np_h20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_np_h20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_np_h20_ca)
summary(model0_tvg_np_h20_ca)
summary.aov(model0_tvg_np_h20_ca)
#Nada

#TX12
tvg_np_h20_tx <- tvg_np_h20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_np_h20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_np_h20_tx)
summary(model0_tvg_np_h20_tx)
summary.aov(model0_tvg_np_h20_tx)
#TxTolyp         7377.3     2469.8   2.987  0.00871 **
#Same overall model NS to begin with 


## NA
tvg_na_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "NA")

model0_tvg_na_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_na_h20)
summary(model0_tvg_na_h20)
summary.aov(model0_tvg_na_h20)
#host.ptgMDCK.TX12(H3N2)  -2738.6     1248.5  -2.193  0.03415 *

#CA09
tvg_na_h20_ca <- tvg_na_h20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_na_h20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_na_h20_ca)
summary(model0_tvg_na_h20_ca)
summary.aov(model0_tvg_na_h20_ca)
#Nada

#TX12
tvg_na_h20_tx <- tvg_na_h20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_na_h20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_na_h20_tx)
summary(model0_tvg_na_h20_tx)
summary.aov(model0_tvg_na_h20_tx)
#TxTolyp       -4163.67    1824.47  -2.282  0.03864 * 
#TxLepto       -5725.33    1824.47  -3.138  0.00726 **

## M
tvg_m_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "M")

model0_tvg_m_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_m_h20)
summary(model0_tvg_m_h20)
summary.aov(model0_tvg_m_h20)
#Nada

## NS
tvg_ns_h20 <-  gcay.screen.count.sgmt %>%
  filter(veh.Tx == "dH2O") %>%
  filter(geno.sgmt == "NS")

model0_tvg_ns_h20 <- lm(count.tvg.sgmt ~ Tx + bioreplicate + host.ptg, data = tvg_ns_h20)
summary(model0_tvg_ns_h20)
summary.aov(model0_tvg_ns_h20)
#host.ptgMDCK.TX12(H3N2) -13377.1     2883.1  -4.640  3.4e-05 ***

#CA09
tvg_ns_h20_ca <- tvg_ns_h20 %>%
  filter(host.ptg == "MDCK.CA09(H1N1)")
model0_tvg_ns_h20_ca <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_ns_h20_ca)
summary(model0_tvg_ns_h20_ca)
summary.aov(model0_tvg_ns_h20_ca)
#Nada

#TX12
tvg_ns_h20_tx <- tvg_ns_h20 %>%
  filter(host.ptg == "MDCK.TX12(H3N2)")
model0_tvg_ns_h20_tx <- lm(count.tvg.sgmt ~ Tx + bioreplicate, data = tvg_ns_h20_tx)
summary(model0_tvg_ns_h20_tx)
summary.aov(model0_tvg_ns_h20_tx)
#TxTolyp         5766.3     2720.8   2.119   0.0500 .

######## Testing Alpelisib on its own ############

alpelisib <- c("dH2O", "Alpe")
#TVG
gcay.screen.count.sgmt %>%
  filter(Tx %in%  alpelisib) %>%
  #filter(host.ptg == "MDCK.CA09(H1N1)")  %>%
  mutate(geno.sgmt = fct_relevel(geno.sgmt, c("PB2", "PB1", "PA",  "HA",  "NP",  "NA",  "M", "NS")))  %>%
  ggplot(aes(x=geno.sgmt, y=count.tvg.sgmt, color=bioreplicate)) +
  #ylim(0, 1.0) +
  xlab("Viral Genome Segment") +
  #ylab("Proportion of Defective Viral Genomes (DelVGs)") +
  theme_tufte() + 
  theme(text = element_text(size = 18,  family="Helvetica")) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(panel.grid.major.y = element_line(color = "lightgray",size = 0.5)) +
  geom_point(aes(fill = bioreplicate), pch=21, size = 3, colour = "black") + theme(legend.position = "none") + facet_grid(host.ptg~Tx, scales = "free")
