#########################################################
#Analysis for publication of Target_Null_Trial
#Sam Forbes
#########################################################

#Load packages
library(lme4)
library(lmerTest)
library(dplyr)
library(xtable)
library(ggplot2);theme_set(theme_classic(base_size = 28))
library(eyetrackingR)

tnt <- readRDS('Data/TNT_All.rds')
tnt$TrialType <- as.character(tnt$TrialType)

tnt <- subset(tnt, TrialType != 'Filler')

tnt$Extra[tnt$TrialType == 'Unrelated'] <- 'Unrelated'
tnt$Extra[tnt$TrialType == 'Knowledge'] <- 'Target'
tnt$Extra[is.na(tnt$Extra)] <- 'Related'

names(tnt)[3] <- 'Trial'

tnt$Other[is.na(tnt$Other)] <- 1
tnt$Other <- as.numeric(tnt$Other)

tnt$Target[is.na(tnt$Target)] <- 0
tnt$Distractor[is.na(tnt$Distractor)] <- 0

workdata <- subset(tnt, Timestamp >= 2000)
workdata <- subset(workdata, TrialType == 'Colour' |
                     TrialType == 'Semantic')

################################################################################
#Original EJ Plot
###############################################################################
tnt2 <- subset(tnt, Timestamp >=2000)

work <- make_eyetrackingr_data(tnt2,
                                participant_column = 'ID',
                                trial_column = 'Trial',
                                time_column = 'Timestamp',
                                trackloss_column = 'Other',
                                aoi_columns = c('Target','Distractor'),
                                treat_non_aoi_looks_as_missing = T)

res <- subset_by_window(work,
                        window_start_time = 2000,
                        window_end_time = 4000,
                        rezero = TRUE,
                        remove = FALSE)

res_time <- make_time_sequence_data(res,
                                    time_bin_size = 100,
                                    predictor_columns = c('Extra'),
                                    aois = 'Target',
                                    summarize_by = 'ID')

EJ <- ggplot(data = res_time,
             aes(x = Time,
                 y = Prop,
                 colour = Extra,
                 linetype = Extra,
                 shape = Extra)) +
  theme(legend.position = c(0.8, .89)) +
  stat_summary(geom = 'pointrange', fun.data = 'mean_se', size = 1.2) +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c('navyblue', 
                                'turquoise',
                                'grey'))  +
  scale_linetype_manual(values = c(1,3,5)) +
  scale_shape_manual(values = c(16,17,18)) +
  ggsave(filename = 'EJ.pdf', path = 'Plots/', 
         width = 297, height = 210, units = "mm")
EJ
################################################################################
#Reshape to match CDI
################################################################################
part <- read.csv("Data/Participants_Expanded_TNT_FP.csv")
part <- subset(part, TNT == 1)
names(part)[1] <- 'ID'

part4 <- select(part, one_of('ID', 'AgeTest', 'Gender', 'Group',
                             'Crocodile', 
                             'Frog', 'Elephant', 'Cheese', 'Banana', 'Carrot', 
                             'Pig', 'Sheep', 'Milk', 'Chocolate', 'Monkey', 'Strawberry'))

library(reshape2)
part5 <- melt(part4, id = 1:4)
names(part5)[5] <- 'ObjectName'

part5$value[is.na(part5$value)] <- 0

merged <- merge(workdata,part5, by = c('ID','ObjectName'))
names(merged)[18] <- 'ObjectUnderstood'

merged$ObjectUnderstood[merged$ObjectUnderstood > 0] <- 1
table(merged$ObjectUnderstood)

#work out how many trials were removed
merged2 <- subset(merged, ObjectUnderstood == 0)
  mergedagg <- aggregate(Timestamp ~ Trial + ID, sum, data = merged2)
 mergedagg <- aggregate(Trial ~ ID, length, data = mergedagg)
sum(mergedagg$Trial)

  
thing <- aggregate(ObjectUnderstood ~ ID + Trial,
                   sum, data = merged)
table(thing$ObjectUnderstood)

accept <- subset(merged, ObjectUnderstood == 1)

accept$Age_C <- ifelse(accept$Age == '24', .5, -.5)
accept$Cond_C <- ifelse(accept$TrialType == 'Colour', .5, -.5)
accept$Type_C <- ifelse(accept$ObjectType == 'Food', .5, -.5)

############################################################################
#CDI Addition
############################################################################

part$cdi <- apply(part[34:449], 1, function(x) length(which(x>=1)))
part$cdiprop <- part$cdi/416
library(dplyr)
library(reshape2)
part2 <- select(part, one_of('ID', 'AgeTest', 'cdiprop', 'Gender', 'Group',
                             'Blue', 'Black', 'Yellow', 'White', 
                             'Green', 'Red', 'Orange', 'Pink',
                             'Brown', 'Grey', 'Purple', 'Aqua'))

everyother <- select(part, c(1, 6, 592, 3, 7, 34:449))
allother <- melt(everyother, id = 1:5)


part3 <- melt(part2, id = 1:4)

part3$value[part3$value > 0] <- 1

names(part3)[5] <- 'Colour'

big <- merge(accept, part3, by = c('ID', 'Colour'))

big$Knowledge[big$value == 1] <- 'Known'
big$Knowledge[big$value == 0] <- 'Unknown'

big$Knowledge <- as.factor(big$Knowledge)

big$Knowledge <- relevel(big$Knowledge, ref = 'Unknown')

big <- big[!is.na(big$Knowledge),]
################################################################################
#Make Eyetracking R
###############################################################################
ETtnt <- make_eyetrackingr_data(big,
                                participant_column = 'ID',
                                trial_column = 'Trial',
                                time_column = 'Timestamp',
                                trackloss_column = 'Other',
                                aoi_columns = c('Target','Distractor'),
                                treat_non_aoi_looks_as_missing = T)

response_tnt <- subset_by_window(ETtnt,
                                 window_start_time = 2000,
                                 window_end_time = 4000,
                                 rezero = TRUE,
                                 remove = FALSE)

trackloss <- trackloss_analysis(data = response_tnt)

trackloss_subjects <- unique(trackloss[,c('ID','TracklossForParticipant')])

#remove files with trackloss
response_clean <- clean_by_trackloss(data = response_tnt,
                                     trial_prop_thresh = 0.6)

#(removed 26 trials, removed 0 participants)
trackloss_clean <- trackloss_analysis(data = response_clean)
trackloss_clean_subjects <- unique(trackloss_clean[, c('ID','TracklossForParticipant')])

#find out trackloss
trackloss_mean<-mean(1 - trackloss_clean_subjects$TracklossForParticipant)
trackloss_sd<-sd(1- trackloss_clean_subjects$TracklossForParticipant)

#See trials contributed by each participant
final_summary <- describe_data(response_clean, 'Target', 'ID')
mean_num_trials<-mean(final_summary$NumTrials)
sd_sum_trials<-sd(final_summary$NumTrials)

#See total attention
trackloss_clean_subjects$ID <- as.character(trackloss_clean_subjects$ID)
trackloss_clean_subjects$Age <- substr(trackloss_clean_subjects$ID, 1, 2)

attention_table <- trackloss_clean_subjects %>%
  group_by(Age) %>%
  summarise(`Mean Attention` = 1-mean(TracklossForParticipant), sd = sd(1-TracklossForParticipant))

xtable(attention_table)
#response_clean$Target[is.na(response_clean$Target)] <- 'FALSE'

response_time <- make_time_sequence_data(response_clean,
                                         time_bin_size = 100,
                                         predictor_columns = c('Age_C',
                                                               'Age',
                                                               'Knowledge',
                                                               'TrialType',
                                                               'TrialType2',
                                                               'cdiprop'),
                                         aois = 'Target',
                                         summarize_by = 'ID')
#write.csv(response_clean,'response_clean2')

response_clean2 <- subset_by_window(response_clean,
                                    window_start_time = 0,
                                    window_end_time = 3000,
                                    remove = T)

response_window <- make_time_window_data(response_clean2,
                                         predictor_columns = c('Age_C',
                                                               'Age',
                                                               'Knowledge',
                                                               'TrialType',
                                                               'TrialType2',
                                                               'cdiprop'),
                                         aois = 'Target')


###############################################################################
#EJ style plot
###############################################################################
EJ2 <- subset(response_time, TrialType2 != 'Semantic') %>%
  ggplot(aes(x = Time,
             y = Prop,
             colour = TrialType2,
             linetype = TrialType2, shape = TrialType2)) +
  stat_summary(geom = 'pointrange', #
               fun.data = 'mean_se', 
               size = 1.2) +
  theme(legend.position = c(0.8, .89)) +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = c('navyblue', 
                                 'turquoise',
                                 'grey'),
                      labels= c('Colour Related',
                                'Colour Only')) +
  scale_linetype_manual(values = c(1,3,5), 
                        labels= c('Colour Related',
                                  'Colour Only')) +
  scale_shape_manual(values = c(16,17),
                     labels= c('Colour Related',
                               'Colour Only')) +
  ggsave(filename = 'EJ_Colour_Only.pdf', path = 'Plots/', 
         width = 297, height = 210, units = "mm")

EJ2
response <- subset(response_time, Time <= 3000)

response2 <- subset(response, SamplesTotal != 0)

colouronly <- subset(response_window, TrialType2 != 'Semantic')




###################################################################
#Plots
####################################################################

bigagg <- aggregate(Prop ~ Age  + ID + TrialType ,
                    mean, data = response2)
bigagg3 <- aggregate(Prop ~ Age  + ID  + TrialType + TrialType2,
                     mean, data = response2)
coll <- ggplot(bigagg,
               aes(x = TrialType,
                   y = Prop,
                   group = TrialType,
                   fill = TrialType,
                   linetype = TrialType)) +
  geom_boxplot(alpha = 0.7) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  scale_fill_manual(values = c('navyblue', 'turquoise'))

coll

t.test(Prop ~ TrialType, data = bigagg)
cc <- subset(bigagg, TrialType == 'Colour')
t.test(cc$Prop, mu = 0.5)
ss <- subset(bigagg, TrialType == 'Semantic')
t.test(ss$Prop, mu = 0.5)

cc4 <- subset(bigagg3, TrialType == 'Colour')
col12 <- ggplot(cc4,
                aes(x = TrialType2,
                    y = Prop,
                    group = TrialType2,
                    fill = TrialType2,
                    linetype = TrialType2)) +
  geom_boxplot(alpha = 0.7) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  scale_fill_manual(values = c('navyblue', 'turquoise'))
col12


t.test(Prop ~ TrialType2,  data = cc4)
#Justify collapsing colour trials

ss$Age2 <- as.character(ss$Age)
sem12 <- ggplot(ss,
                aes(x = Age2,
                    y = Prop,
                    group = Age2)) +
  geom_boxplot(alpha = 0.7) +
  theme_classic(base_size = 30) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  xlab('Participant Age') +
  ggsave(filename = 'Semantic_Only.pdf', path = 'Plots/', 
         width = 297, height = 210, units = "mm")
sem12

t.test(Prop ~ Age2, data = ss)
#################################################################################
#Modelling bu age
################################################################################

#Load packages
library(lme4)
library(lmerTest)
response2$Knowledge_C <- ifelse(response2$Knowledge == 'Known', .5, -.5)
bigcol <- subset(response2, TrialType == 'Colour')



agemodfull <- glmer(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
                      (ot1 + ot2 + ot3 + ot4) *
                      Age_C +
                      (1 + ot1||ID),
                    family = binomial,
                    na.action = na.exclude,
                    glmerControl(optimizer = 'bobyqa'), 
                    #optCtrl = list(maxfun = 20000)),
                    data = bigcol)
save(agemodfull, file = 'Data/agemodfull')

bigcol$colpred <- predict(agemodfull, type = 'response')
bigcol$Age <- as.character(bigcol$Age)
bigcol$Age <- as.factor(bigcol$Age)
bigcol$Age <- relevel(bigcol$Age, '24')

a <- ggplot(bigcol,
            aes(x = Time,
                y = Prop,
                linetype = Age,
                shape = Age,
                colour = Age)) +
  stat_summary(aes(y = colpred),
               geom = 'line', #
               fun.y = 'mean', 
               size = 1) +
  stat_summary(geom = 'pointrange',
               fun.data = 'mean_se',
               size = 1)+
  theme_classic(base_size = 24) +
  theme(legend.position = c(0.3, .89),
        legend.key.size = unit(1.4, 'lines')) +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  scale_colour_manual(values = c('navyblue', 
                                 'turquoise')) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_shape_manual(values = c(16,17,18))


#################################################################################
#By Knowledge
#################################################################################
knowmodfull <- glmer(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
                       (ot1 + ot2 + ot3 + ot4) *
                       Knowledge_C +
                       (1 + ot1 ||ID),
                     family = binomial,
                     na.action = na.exclude,
                     glmerControl(optimizer = 'bobyqa'), 
                     #optCtrl = list(maxfun = 20000)),
                     data = bigcol)
save(knowmodfull, file = 'Data/knowmodfull')

anova(agemodfull, knowmodfull, test = 'chisq') #for AIC
#           Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
#agemodfull  12 85375 85457 -42676    85351                             
#knowmodfull 12 85294 85375 -42635    85270 81.457      0  < 2.2e-16 ***

bigcol$knowpred <- predict(knowmodfull, type = 'response')
bigcol$Knowledge <- as.character(bigcol$Knowledge)

k <- ggplot(bigcol,
            aes(x = Time,
                y = Prop,
                linetype = Knowledge,
                shape = Knowledge,
                colour = Knowledge)) +
  stat_summary(aes(y = knowpred),
               geom = 'line', #
               fun.y = 'mean', 
               size = 1) +
  stat_summary(geom = 'pointrange',
               fun.data = 'mean_se',
               size = 1)+
  theme_classic(base_size = 24) +
  theme(legend.position = c(0.3, .89),
        legend.key.size = unit(1.4, 'lines')) +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  scale_colour_manual(name = 'Colour Word Knowledge',
                      values = c('navyblue', 
                                 'turquoise')) +
  scale_shape_manual(name = 'Colour Word Knowledge',
                     values = c(16, 17)) +
  scale_linetype_manual(name = 'Colour Word Knowledge',
                        values = c(1, 2))

source('Misc/multiplot.R')
multiplot(a,k) +
ggsave(filename = 'AgeKnowledge.pdf', path = 'Plots/', 
       width = 297, height = 210, units = "mm")

nineteen <- subset(bigcol, Age == '19')
################################################################################
#nineteen months
################################################################################

#nine <- read.csv('nine.csv')
#nine[1] <- NULL



ninemodfull <- glmer(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
                       (ot1 + ot2 + ot3) *
                       Knowledge_C +
                       (1 + ot1 ||ID),
                     family = binomial('logit'),
                     na.action = na.exclude,
                     glmerControl(optimizer = 'bobyqa', 
                                  optCtrl = list(maxfun = 20000)),
                     data = nineteen)
save(ninemodfull, file = 'Data/ninemodfull')

nineteen$Predflat <- predict(ninemodfull, type = 'response')
nineteen$Knowledge <- as.character(nineteen$Knowledge)

ggplot(nineteen,
       aes(x = Time,
           y = Prop,
           linetype = Knowledge,
           shape = Knowledge,
           colour = Knowledge)) +
  stat_summary(aes(y = Predflat),
               geom = 'line', #
               fun.y = 'mean', 
               size = 1) +
  stat_summary(geom = 'pointrange',
               fun.data = 'mean_se',
               size = 1.2)+
  theme(legend.position = c(0.4, .85)) +
  geom_hline(yintercept = .5, linetype = 'dashed',colour = 'black') +
  ylab('Proportion Looking') +
  scale_colour_manual(name = 'Knowledge of \nColour Word',
                      values = c('navyblue', 
                                 'turquoise')) +
  scale_shape_manual(name = 'Knowledge of \nColour Word',
                     values = c(16, 17)) +
  scale_linetype_manual(name = 'Knowledge of \nColour Word',
                        values = c('solid','longdash')) +
ggsave(filename = 'Nineteen.pdf', path = 'Plots/', 
       width = 297, height = 210, units = "mm")

######################################################################
#Chi-Squared
######################################################################
big3 <- aggregate(Prop ~ ID  + Age + Knowledge, mean, data = bigcol)
tbl <- table(big3$Age, big3$Knowledge)
chisq.test(tbl)

######################################################################
#Other 19 m.o. comparisons
######################################################################
#only include those with cdi answers
nine2 <- subset(nineteen, cdiprop>0)


ninemodcomp <- glmer(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
                       (ot1 + ot2 + ot3) *
                       Knowledge_C +
                       (1 + ot1 ||ID),
                     family = binomial('logit'),
                     na.action = na.exclude,
                     glmerControl(optimizer = 'bobyqa', 
                                  optCtrl = list(maxfun = 20000)),
                     data = nine2)
save(ninemodcomp, file = 'Data/ninemodcomp')

ninemodcdi <- glmer(cbind(SamplesInAOI, SamplesTotal-SamplesInAOI) ~
                      (ot1 + ot2 + ot3) *
                      cdiprop +
                      (1 + ot1 ||ID),
                    family = binomial('logit'),
                    na.action = na.exclude,
                    glmerControl(optimizer = 'bobyqa', 
                                 optCtrl = list(maxfun = 20000)),
                    data = nine2)
save(ninemodcdi, file = 'Data/ninemodcdi')

#Obviously cannot be compared with chisq, but just for AIC
anova(ninemodcdi,ninemodcomp)

###################################################################
#All other words (SIM)
###################################################################
other <- subset(everyother, cdiprop > 0)
other <- subset(other, Group == '19')

other[is.na(other)] <- 0
other[c(6:421)] <- ifelse(other[c(6:421)] > 0, 1, 0)

nine3 <- subset(nine2, Time <= 2500 & Time >= 1500)
nine3 <- aggregate(cbind(SamplesInAOI, SamplesTotal) ~ ID +
                     Knowledge + Knowledge_C, sum, data = nine3)
nine3$Prop <- nine3$SamplesInAOI / nine3$SamplesTotal

sims <- merge(nine3, other, by = 'ID')

#convert to character
for(i in c(11:426)) {
  sims[,i] <- as.character(sims[,i])
}

# Set participants as looking or not
sims$Looked[sims$Prop > 0.55] <- 'Looked'
sims$Looked[sims$Prop <= 0.55] <- 'No Look'
sims$Looked <- as.character(sims$Looked)

tab <- table(sims$Knowledge, sims$Looked) 

#remove nans
sims <- sims[!is.nan(sims$Prop),]

#empty data frame
props <- rep(0, 416)

#run a loop
for (i in c(11:426)) {
  vari <- as.numeric(sims[[i]])
  y <- subset(sims, sims[i] == '1')
  a <- mean(y$Prop)
  j <- i - 10
  props[j] <- a
} 

simsk <- subset(sims, Knowledge == 'Known')
mn <- mean(simsk$Prop)

t.test(props, mu = mn )

props <- data.frame(props)

maincols <- c(337, 350, 362, 371)

props$Word <- NA
props$Word[maincols] <- c('blue', 'green', 'red', 'yellow')
props$height <- NA
props$height[maincols] <- c(30, 45, 60, 75)

library(ggrepel)
ggplot(props,
       aes(x = props, label = Word)) +
  geom_histogram(fill = 'blue',
                 alpha = 0.4, colour = NA,
                 binwidth = 0.02) +
  geom_vline(xintercept = mn,
             colour = 'red',
             linetype = 'dashed',
             size = 1) +
  geom_text_repel(aes(y = 30), nudge_y = 5, force = 100,
                  min.segment.length = 0, angle = 90,  size = 7) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, .89)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,130)) +
  ylab('Count') +
  xlab('Prop looks to target') +
  ggsave(filename = 'Simulated.pdf', path = 'Plots/', 
         width = 297, height = 210, units = "mm")




