library(readxl)
library(ggplot2)
library(survival)
library(ggpubr)
library(survminer) 
library(bdsmatrix) 
library(coxme)
library(plyr) 
library(lattice)
library(Rmisc)
library(emmeans)
library(dplyr) 
library(carData)
library(car) 
library(ARTool)
library(emmeans)
library(tidyr)
library(forecast)
library(zoo) 
library(gridExtra)
setwd("~/Library/CloudStorage/OneDrive-UCB-O365/Umich/Umich R files/Thesis")

###################### SET UP ########################################
 setwd("~/Library/CloudStorage/OneDrive-UCB-O365/Umich/Umich R files/Thesis")
rm(list=ls())
###################### GUT PENETRATION ########################################
gut_penetration <- read.csv("gut_penetration_mkd.csv") #load data
gut_penetration <- subset(gut_penetration, exptime== "early")
gut_penetration$spores <- (gut_penetration$embup + gut_penetration$embdw + gut_penetration$penup + gut_penetration$pendw)
gut_penetration$spores<-as.numeric(gut_penetration$spores) # change the "spores" data into values

hist(gut_penetration$spores,
     main = "Histogram of Spore Counts",
     xlab = "Spores",
     ylab = "Frequency",)

M1<-glm(spores~parasite*food, data = gut_penetration) # creates linear model 
anova(M1) # prints stats



average_attacking <- gut_penetration %>%
  group_by(food) %>%
  summarize(attaching = mean(spores, na.rm = TRUE))
print(average_attacking)

op <- par(mfrow = c(2, 2)) # following three lines create diagnostic plots on non square rooted value
plot(M1)             
par(op) 


gptemp<-na.omit(gut_penetration) # removes all NAs from the database (necessary for summarySE)
gptempSE<-summarySE(data=gptemp,measurevar = 'spores', groupvars = c('food','parasite'))



SVP2 <- ggplot(gut_penetration, aes(x = food, y = spores, fill = interaction(parasite, food))) + 
    geom_violin(position = position_dodge(.8), alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
                 position = position_dodge(.8), color = "white") +
    scale_fill_manual(values = c('#09753f','red','blue','orange'),
                      labels = c("AM High","OP/AM High","AM Low","OP/AM Low"),
                      name = "Pathogen") +
    ylab("Attacking spores") +
    xlab("Food") +
    scale_y_continuous(limits = c(0, 32), expand = c(0, 0)) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)
    )

        
ggsave("SVP2_MKD.png",SVP2, width = 8, height = 8) # exports the plot to .png
SVP2

#spore_plots <- grid.arrange(SVP2$plot$SVP2, SVP1$plot$SVP1, ncol = 2) #arrange the plots

#ggsave("spore_plots.png", plot = spore_plots, width = 8, height = 4) # exports the plot

################# MORTALITY ####################################################

mortality <- read.csv("mortality_mkd.csv") #load data
mortality$lifespan<-as.numeric(mortality$lifespan)
#mortality$metsch_infection<-as.factor(mortality$metsch_infection)
mortality$food<-as.factor(mortality$food)
mortality$parasite<-as.factor(mortality$parasite)

#mortality_subset_inf <- subset(mortality, 
 #                              parasite %in% c("Metsch", "MicG/Metsch") & metsch_infection == 1)
#mortality_subset_uninf <- subset(mortality, 
 #                              parasite %in% c("Metsch", "MicG/Metsch") & metsch_infection == 0)

mortality_filtered <- subset(mortality, 
                                 parasite %in% c("MET", "MG_MET"))
                                 
mortality_filtered$parasite <- droplevels(mortality_filtered$parasite)

average_lifespan <- mortality_filtered %>%
  group_by(food, parasite, ) %>%
  summarize(average_lifespan = mean(lifespan, na.rm = TRUE))
print(average_lifespan)

M_lifespan <- glm(lifespan ~ parasite * food, data = mortality_filtered, family = poisson) # creates linear model 
anova_results <- anova(M_lifespan, test = "Chisq")
print(anova_results)
summary(M_lifespan)

op <- par(mfrow = c(2, 2)) # following three lines create diagnostic plots 
plot(M_lifespan)             
par(op) 

hist(mortality_filtered$lifespan,
     main = "Histogram of infected mortality",
     xlab = "lifespan",
     ylab = "Frequency",)

mortality_filtered<-na.omit(mortality_filtered)      # removing the rows with NA

cox<- coxph(Surv(lifespan) ~ parasite*food, data = mortality_filtered) # setting up cox model
Anova(cox) # prints stats



LM_filtered_MKD<- # produces Kaplan-Meier plot for Ank Food treatments
  ggsurvplot(survfit(Surv(lifespan,   status) ~ parasite +food, data = mortality_filtered),
             conf.int=F,
             xlim = c(0, 65)
             )
LM_filtered_MKD$plot<-LM_filtered_MKD$plot+
  ggplot2::annotate("text", 
                    x = 10, y = 0.2, # x and y coordinates of the text
                    label = "",size = 10)+
  scale_color_manual(values = c('#09753f','blue','red','orange')
                     , 
                     labels=c("AM High","AM Low","OP/AM High","OP/AM Low")
                     )+ 
  xlab("Day")+
  ylab("")+
  ggtitle("Experiment 1\n")+
  scale_y_continuous(limits=c(0,1.05),expand = c(0,0))+
  scale_x_continuous(limits = c(0,90), expand = c(0,0))+
  theme(legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.7,0.5),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(color = "black",size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15))+
    geom_text(aes(x = 4, y = 0.9, label = "A)"), size = 6)

ggsave("Kaplan_Meier_Plot_MKD.png", plot = LM_filtered_MKD$plot, width = 8.5, height = 6, units = "in", dpi = 300)
LM_filtered_MKD


############################ REPRODUCTION ######################################

reproduction <- read.csv("reproduction_mkd.csv") #load data
#totbabies_metsch <- read.csv("totbabies_metsch.csv") #load data
reproduction$age<-as.numeric(reproduction$age) # turning age and NB into numeric
reproduction$NB<-as.numeric(reproduction$NB)
reproduction$rep<-as.character(reproduction$rep)
reproduction$parasite<-as.factor(reproduction$parasite)
reproduction$food<-as.factor(reproduction$food)
#totbabies_metsch$metsch_infection <- as.factor(totbabies_metsch$metsch_infection)
#totbabies_metsch$parasite <- factor(totbabies_metsch$parasite)
#totbabies_metsch$food <- factor(totbabies_metsch$food)
#reproduction<-na.omit(reproduction) # removes units with NA

str(reproduction)
summary(reproduction)


reproduction <- subset(reproduction, parasite == "MET" | parasite == "MG_MET")

hist(reproduction$NB,
     main = "Histogram of totbabies",
     xlab = "babies",
     ylab = "Frequency",)

totbabies <- reproduction %>% group_by(food,parasite,rep) %>% # next four lines sum all the offspring
  summarise(totbabies = sum(NB),                                       # and gives total lifetime reproduction
          .groups = 'drop'
            )                                       # as totbabies variable in new dataframe
totbabies <- totbabies %>% as.data.frame()                        # named totbabies
totbabies<-na.omit(totbabies) # removes units with NA


average_totbabies <- totbabies %>%
  group_by(food, parasite) %>%
  summarize(average_totbabies = mean(totbabies, na.rm = TRUE))
print(average_totbabies)

hist(totbabies$totbabies, breaks=30)

repM <- glm(totbabies ~ parasite * food, family = poisson, data = totbabies)
op <- par(mfrow = c(2, 2))                         # diagnostic plots
plot(repM)             
par(op)
Anova(repM, type = 2)

totbabies_SE<-summarySE(data = totbabies, measurevar = "totbabies", groupvars = c("food","parasite"), na.rm=T)

# Create the plot with the filtered data
R3 <- ggplot(totbabies_SE, aes(
  x = food, y = totbabies,
  fill = interaction(parasite, food),
  shape = interaction(parasite, food)
)) + 
  geom_errorbar(aes(ymin=totbabies-se, ymax=totbabies+se), colour="black", width=.3, position=position_dodge(0.4)) +
  geom_point(position=position_dodge(0.4), size=4)+
  scale_y_continuous(limits = c(0,40), expand = c(0,0))+
  xlab("Food")+
  ylab("Lifetime offspring production")+
  ggtitle("Experiment 1\n")+
scale_x_discrete(labels=c("High","Low"))+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(size = 5, colour = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.position = c(0.8, 0.8) 
  ) +
  scale_fill_manual(
    name = "",
    labels = c("AM High", "OP/AM High", "AM Low", "OP/AM Low"),
    values = c("#09753f", "red", "blue", "orange")
  ) +
  scale_shape_manual(
    name = "",
    labels = c("AM High", "OP/AM High", "AM Low", "OP/AM Low"),
    values = c(21, 24, 21, 24)
  ) +
  geom_text(aes(x = 0.6, y = 38, label = "A)"), size = 6) 

ggsave("Reproduction_MKD.png", plot = R3, width = 4, height = 4) # exports the plot
R3


################################################################################

####################### PREVALENCE AND BURDEN ##################################

metsch_infection <- read.csv("metsch_infection_MKD.csv") #load data
#metsch_infection$Mature<-as.numeric(metsch_infection$mature) # changes n. of mature spores into numeric
#metsch_infection$infected<-as.numeric(metsch_infection$infected) # changes effective infection score into numeric (not sure what this is)
#metsch_infection$MetschAll<-as.numeric(metsch_infection$MetschAll) # changes infection score into numeric
metsch_infection$food<-as.factor(metsch_infection$food) # changes food into factor
metsch_infection$parasite<-as.factor(metsch_infection$parasite) # changes parasite into factor
#metsch_infection<-na.omit(metsch_infection) # removes all NA cases

#### INFECTION PROBABILITY ####
subset_metsch_exposure <- subset(metsch_infection,(parasite == "MET" | parasite == "MG_MET"))

prevM1<-glm(MetschAll~parasite*food,family = binomial(link = "logit"),data = subset_metsch_exposure) # global binomial glm model with interaction
prevM1
Anova(prevM1)


prevSE<-summarySE(data = subset_metsch_exposure, measurevar = "MetschAll", groupvars = c("parasite", "food"), na.rm = T)
# the line above summarizes data within treatments and calculates sd, se and ci

S2<- # makes the plot
  ggplot(prevSE, aes(x=food, y=MetschAll, fill=interaction(parasite,food), shape = interaction(parasite,food),
                     #alpha = food
                     )) + 
  geom_line(position=position_dodge(.4))+
  scale_alpha_discrete(range = c(1, 0.5))+
  geom_errorbar(aes(ymin=MetschAll-ci, ymax=MetschAll+ci), width=.2, 
                position=position_dodge(.4)) + 
  scale_y_continuous(limits = c(-0.041,1.1))+
  ylab("A. monospora \n infection probability")+
  xlab("Food")+
  ggtitle("Experiment 1\nInfection Probability")+
  scale_x_discrete(labels=c("High","Low"))+
  geom_point(position=position_dodge(.4), size=4)+
  theme_classic()+
  theme_classic() +
  theme(plot.margin = margin(20, 0, 0, 0),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 15, colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.title = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 10, color = "black"),
    legend.position = c(0.7,0.3) 
  ) +
    scale_fill_manual(
        name = "",
        labels = c("AM High", "OP/AM High", "AM Low", "OP/AM Low"),
        values = c("#09753f", "red", "blue", "orange")
    ) +
    scale_shape_manual(
        name = "",
        labels = c("AM High", "OP/AM High", "AM Low", "OP/AM Low"),
        values = c(21, 24, 21, 24) 
    ) +
    geom_text(aes(x = 0.6, y = 0.95, label = "A)"), size = 6)
S2 # prints the plot

ggsave("metsch_infection_prob_MKD.png", plot = S2, width = 4, height = 4) # exports the plot

#### MATURE SPORES BURDEN ####
InfNoZero <- subset(subset_metsch_exposure, !is.na(MetschAll) & M > 0)
#InfNoZero<-na.omit(InfNoZero) # removes cases with NA
InfNoZero$parasite <- factor(InfNoZero$parasite)
InfNoZero$food <- factor(InfNoZero$food)
InfNoZero$M <- as.numeric(InfNoZero$M)
#InfNoZero<-InfNoZero[,-which(names(InfNoZero) == "X")]
#InfNoZero<-InfNoZero[,-which(names(InfNoZero) == "X.1")]

par(mfrow=c(1,1))
hist(InfNoZero$M,
     main = "Histogram of Metsch Spores",
     xlab = "Spores",
     ylab = "Frequency",
     breaks=10)

mature_burden2<-lm(M~parasite*food,data=InfNoZero) # creates a linear model
summary(mature_burden2) # prints stats


op <- par(mfrow = c(2, 2)) 
plot(mature_burden2)             
par(op)


metsch_infection_matureSE<-summarySE(data = InfNoZero, measurevar = "M", groupvars = c("parasite", "food"), na.rm=T)
# the line above summarizes data within treatments and calculates sd, se and ci

# Create the plot with all 4 groups for consistency but the analysis is on treatment
S3 <- ggplot(metsch_infection_matureSE, aes(
  x = food, y = M,
  fill = interaction(parasite, food),
  shape = interaction(parasite, food)
)) + 
  geom_errorbar(aes(ymin=M-se, ymax=M+se), colour="black", width=.3, position=position_dodge(0.4)) +
  geom_point(position=position_dodge(0.4), size=4)+
  scale_y_continuous(limits = c(0,100), expand = c(0,0))+
  ylab("A. monospora mature spore\nburden x10,000")+
  xlab("Food")+
  ggtitle("Experiment 1\nSpore Burden")+
    scale_x_discrete(labels=c("High","Low"))+
  theme_classic()+
  theme(plot.margin = margin(20, 0, 0, 0),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.title = element_text(size = 5, colour = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.position = c(0.3,0.3) 
  ) +
  scale_fill_manual(
    name = "",
    labels = c("AM High", "OP/AM High", "AM Low", "OP/AM Low"),
    values = c("#09753f", "red", "blue", "orange")
  ) +
  scale_shape_manual(
    name = "",
    labels = c("AM High", "OP/AM High", "AM Low", "OP/AM Low"),
    values = c(21, 24, 21, 24)
  )  +
    geom_text(aes(x = 0.6, y = 92, label = "C)"), size = 6)

ggsave("spores_mkd.png", plot = S3, width = 4, height = 4) # exports the plot
S3

#ggsave("spore_plots_mkd.png", plot = spore_plots, width = 8, height = 4) # exports the plot
