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

###################### GUT PENETRATION ########################################
gut_penetration <- read.csv("penetration.csv") #read-in data

hist(gut_penetration$spores_piercing, #visualize distribution of data in a histogram
     main = "Histogram of Spore Counts",
     xlab = "Spores",
     ylab = "Frequency",)

M2<-glm(spores_piercing~parasite*food, data = gut_penetration, family= poisson) # creates generalized linear model 
summary(M2)# prints stats


average_attacking <- gut_penetration %>% #calcualtes average number of attacking spores by group
  group_by(food) %>%
  summarize(attaching = mean(spores_piercing, na.rm = TRUE))
print(average_attacking)

op <- par(mfrow = c(2, 2)) # creates diagnostic plots
plot(M2)             
par(op) 


gptemp<-na.omit(gut_penetration) # removes all NAs from the database (necessary for summarySE)
gptempSE<-summarySE(data=gptemp,measurevar = 'spores_piercing', groupvars = c('food','parasite')) #summary statistics for attacking spores data


SVP1 <- ggplot(gut_penetration, aes(x = food, y = spores_piercing, fill = interaction(parasite, food))) + #creates violin plot for the attacking spores data
  geom_violin(position = position_dodge(.8), alpha = 0.5) +
  stat_summary(fun=mean, geom="point", shape=23, size=2, position=position_dodge(.8), color = "white")+
  scale_fill_manual(values = c('#09753f','red','blue','orange'), labels=c("MB Ank","OP/MB Ank","MB CYA","OP/MB CYA"), name = "Pathogen")+
  ylab("Attacking spores")+
  xlab("Food")+
  theme_classic()+
  scale_y_continuous(limits= c(0,5.05),expand = c(0,0))+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.position="none") 
  
ggsave("SVP1.png",SVP1, width = 8, height = 8) # exports the plot to .png
SVP1 #prints violin plot

################# MORTALITY ####################################################
mortality <- read.csv("mortality.csv") #load data
mortality$lifespan<-as.numeric(mortality$lifespan) #lifespan as numeric
mortality$metsch_infection<-as.factor(mortality$metsch_infection) #infection as factor
mortality$food<-as.factor(mortality$food) # food as factor
mortality$parasite<-as.factor(mortality$parasite) # parasite as factor

mortality_subset_inf <- subset(mortality, #subsets data to include only those exposed to metsch and infected
                               parasite %in% c("Metsch", "MicG/Metsch") & metsch_infection == 1)

mortality_subset_uninf <- subset(mortality, #subsets data to include only those exposed to metsch and uninfected
                               parasite %in% c("Metsch", "MicG/Metsch") & metsch_infection == 0)


average_lifespan <- mortality_subset_inf %>% #calculates average lifespan by parasite and food for infected animals
  group_by(food, parasite, ) %>%
  summarize(average_lifespan = mean(lifespan, na.rm = TRUE))
print(average_lifespan)

average_lifespan <- mortality_subset_inf %>% #calculates average lifespan by parasite and food for infected animals
  group_by(food) %>%
  summarize(average_lifespan = mean(lifespan, na.rm = TRUE))
print(average_lifespan)

mortality_subset_inf$parasite <- droplevels(mortality_subset_inf$parasite) #drops other parasite groups from infected data
mortality_subset_uninf$parasite <- droplevels(mortality_subset_uninf$parasite) #drops other parasite groups from uninfected data

op <- par(mfrow = c(2, 2)) # following three lines create diagnostic plots for infected subset
plot(M_lifespan_inf)             
par(op) 

hist(mortality_subset_inf$lifespan, #visualizes infected subset with a histogram
     main = "Histogram of infected mortality",
     xlab = "lifespan",
     ylab = "Frequency",)

op <- par(mfrow = c(2, 2)) # following three lines create diagnostic plots for uninfected subset 
plot(M_lifespan_uninf)             
par(op) 

hist(mortality_subset_uninf$lifespan, #visualizes uninfected subset with a histogram
     main = "Histogram of uninfected mortality",
     xlab = "lifespan",
     ylab = "Frequency",)

#mortality<-na.omit(mortality)      # removing the rows with NA

cox_inf<- coxph(Surv(lifespan) ~ parasite*food, data = mortality_subset_inf) # setting up cox model
Anova(cox_inf) # prints stats

cox_uninf<- coxph(Surv(lifespan) ~ parasite*food, data = mortality_subset_uninf) # setting up cox model
Anova(cox_uninf) # prints stats

LM_inf <- # produces Kaplan-Meier plot for infected subset
    ggsurvplot( 
  survfit(Surv(lifespan, status) ~ parasite + food, data = mortality_subset_inf),
  conf.int = FALSE,
  xlim = c(0, 25)  # Ensure consistency
)

LM_inf$plot <- LM_inf$plot +
  ggplot2::annotate("text", x = 3, y = 0.9, label = "B)", size = 6) +  # Fix text annotation
  scale_color_manual(values = c('#09753f', 'blue', 'red', 'orange'), 
                     labels = c("AM Ank", "AM CYA", "OP/AM Ank", "OP/AM CYA")) +
  xlab("Day") +
  ggtitle("Experiment 2\nInfected with AM") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +  
  scale_x_continuous(limits = c(0, 25), expand = c(0.05, 0)) +  
  theme(legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.3, 0.21),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15))

# Save the plot
ggsave("Kaplan_Meier_Plot_Inf.png", plot = LM_inf$plot, width = 8, height = 6, units = "in", dpi = 300)
LM_inf

LM_uninf<- # produces Kaplan-Meier plot for uninfected subset
    ggsurvplot(survfit(Surv(lifespan,   status) ~ parasite +food, data = mortality_subset_uninf),
               conf.int=F,
               xlim = c(0, 90)
    )
LM_uninf$plot<-LM_uninf$plot+
    ggplot2::annotate("text", 
                      x = 10, y = 0.2, # x and y coordinates of the text
                      label = "",size = 10)+
    scale_color_manual(values = c('#09753f', 'blue', 'red', 'orange'), 
                       labels = c("AM Ank", "AM CYA", "OP/AM Ank", "OP/AM CYA")) +
   
    xlab("Day")+
    ylab("")+
    ggtitle("Experiment 2\nUninfected with AM")+
    scale_y_continuous(limits=c(0,1.05),expand = c(0,0))+
    scale_x_continuous(limits = c(0,90), expand = c(0,0))+
    theme(legend.text = element_text(size = 10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15),
          legend.position = c(0.21,0.21),
          axis.text.x = element_text(colour = "black", size = 15),
          axis.text.y = element_text(color = "black",size = 15),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(colour = "black", size = 15))+
    geom_text(aes(x = 4, y = 0.9, label = "C)"), size = 6)

ggsave("Kaplan_Meier_Plot_Uninf.png", plot = LM_uninf$plot, width = 8.5, height = 6, units = "in", dpi = 300)
LM_uninf

mortality_plots <- grid.arrange(LM_filtered_MKD$plot, LM_inf$plot, LM_uninf$plot, ncol = 3) #arrange the plots

ggsave("Mortality_plots.png", plot = mortality_plots, width = 12, height = 4) # exports the plot

############################ REPRODUCTION ######################################
rm(list=ls())
totbabies_metsch <- read.csv("totbabies_metsch.csv") #load data
totbabies_metsch$metsch_infection <- as.factor(totbabies_metsch$metsch_infection) #mestch infection as factor
totbabies_metsch$parasite <- factor(totbabies_metsch$parasite) #parasite as factor
totbabies_metsch$food <- factor(totbabies_metsch$food) #food as factor

subset_totbabies_inf <- subset(totbabies_metsch, metsch_infection == "1" & (parasite == "Metsch" | parasite == "MicG/Metsch")) #subsets data to include only those exposed to metsch and infected
subset_totbabies_uninf <- subset(totbabies_metsch, metsch_infection == "0" & (parasite == "Metsch" | parasite == "MicG/Metsch")) #subsets data to include only those exposed to metsch and uninfected

hist(subset_totbabies_inf$totbabies, #visualizes infected subset data with a histogram
     main = "Histogram of totbabies inf",
     xlab = "babies",
     ylab = "Frequency",)

hist(subset_totbabies_uninf$totbabies, #visualizes uninfected data with a histogram
     main = "Histogram of totbabies Uninf",
     xlab = "babies",
     ylab = "Frequency",)



repM_inf <- glm(totbabies ~ parasite * food, family = poisson, data = subset_totbabies_inf) #makes GLM with poisson distribution of infected subset

op <- par(mfrow = c(2, 2)) # diagnostic plots
plot(repM_inf)             
par(op)

Anova(repM_inf, type = 2) #runs Anova on the glm
summary(repM_inf)
?Anova


repM_uninf <- glm(totbabies ~ parasite * food, family = Gamma, data = subset_totbabies_uninf) #makes GLM with poisson distribution of infected subset
op <- par(mfrow = c(2, 2))                        # diagnostic plots
plot(repM_uninf)             
par(op)
 Anova(repM_uninf, type = 2)
emmeans_results_uninf <- emmeans(repM_uninf, pairwise ~ parasite * food) #pairwise comparisons
emmeans_results_uninf #prints pairwise comparisons
summary(emmeans_results_uninf$contrasts, adjust = "tukey")

totbabies_SE_inf<-summarySE(data = subset_totbabies_inf, measurevar = "totbabies", groupvars = c("food","parasite"), na.rm=T)

# creates plot with infected subset
R5 <- ggplot(totbabies_SE_inf, aes(
    x = food, y = totbabies,
    fill = interaction(parasite, food),
    shape = interaction(parasite, food)
)) + 
    geom_errorbar(aes(ymin=totbabies-se, ymax=totbabies+se), colour="black", width=.3, position=position_dodge(0.4)) +
    geom_point(position=position_dodge(0.4), size=4)+
    scale_y_continuous(limits = c(0,40), expand = c(0,0))+
    ylab("")+
    xlab("Food")+
    ggtitle("Experiment 2\nInfected with AM")+
    theme_classic()+
    theme(axis.text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(size = 20, colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 18),
          legend.title = element_text(size = 5, colour = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.position = c(0.3,0.3) 
    ) +
    scale_fill_manual(
        name = "",
        labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
        values = c("#09753f", "red", "blue", "orange")
    ) +
    scale_shape_manual(
        name = "",
        labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
        values = c(21, 24, 21, 24)
    ) +
    geom_text(aes(x = 0.6, y = 38, label = "B)"), size = 6)

ggsave("Reproduction_inf.png", plot = R5, width = 4, height = 4) # exports the plot
R5

totbabies_SE_uninf<-summarySE(data = subset_totbabies_uninf, measurevar = "totbabies", groupvars = c("food","parasite"), na.rm=T)

R4 # creates plot with uninfected subset
R4 <- ggplot(totbabies_SE_uninf, aes(
  x = food, y = totbabies,
  fill = interaction(parasite, food),
  shape = interaction(parasite, food)
)) + 
  geom_errorbar(aes(ymin=totbabies-se, ymax=totbabies+se), colour="black", width=.3, position=position_dodge(0.4)) +
  geom_point(position=position_dodge(0.4), size=4)+
  scale_y_continuous(limits = c(0,350), expand = c(0,0))+
  xlab("Food")+
  ylab("")+
  ggtitle("Experiment 2\nUninfected with AM")+
  theme_classic()+
  theme(axis.text = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_text(size = 5, colour = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.position = c(0.8,0.8) 
  ) +
  scale_fill_manual(
    name = "",
    labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
    values = c("#09753f", "red", "blue", "orange")
  ) +
  scale_shape_manual(
    name = "",
    labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
    values = c(21, 24, 21, 24)
  ) +
  geom_text(aes(x = 0.9, y = 330, label = "a"), size = 4) +
  geom_text(aes(x = 1.1, y = 180, label = "b"), size = 4) +
  geom_text(aes(x = 1.9, y = 170, label = "b"), size = 4) +
  geom_text(aes(x = 2.1, y = 80, label = "c"), size = 4) +
  geom_text(aes(x = 0.6, y = 320, label = "C)"), size = 6) 

ggsave("Reproduction_uninf.png", plot = R4, width = 5, height = 4) # exports the plot
R4



reproduction_plots <- grid.arrange(R3, R5, R4, ncol = 3) #arranges plots
ggsave("Reproduction_plots.png", plot = reproduction_plots, width = 12, height = 4) # exports the plots

####################### PREVALENCE AND BURDEN ##################################

metsch_infection <- read.csv("metsch_infection.csv") #load data
metsch_infection$Mature<-as.numeric(metsch_infection$mature) # changes n. of mature spores into numeric
metsch_infection$infected<-as.numeric(metsch_infection$infected) # changes effective infection score into numeric
metsch_infection$food<-as.factor(metsch_infection$food) # changes food into factor
metsch_infection$parasite<-as.factor(metsch_infection$parasite) # changes parasite into factor
#metsch_infection<-na.omit(metsch_infection) # removes all NA cases

#### INFECTION PROBABILITY ####
subset_metsch_exposure <- subset(metsch_infection,(parasite == "Metsch" | parasite == "MicG/Metsch")) #subsets data into only those exposed to metsch

prevM1<-glm(infected~parasite*food,family = binomial(link = "logit"),data = subset_metsch_exposure) # makes global binomial glm model with interaction
prevM1
Anova(prevM1)

prevSE<-summarySE(data = subset_metsch_exposure, measurevar = "infected", groupvars = c("parasite", "food"), na.rm = T) #summarizes data within treatments and calculates sd, se and ci

S4<- # makes the plot
  ggplot(prevSE, aes(x=food, y=infected, fill=interaction(parasite,food), shape = interaction(parasite,food),
                     #alpha = food
                     )) + 
  geom_line(position=position_dodge(.4))+
  scale_alpha_discrete(range = c(1, 0.5))+
  #guides(fill=guide_legend("Parasite"), alpha = "none")+
  #expand_limits(y=c(0,4))+
  geom_errorbar(aes(ymin=infected-ci, ymax=infected+ci), width=.2, 
                position=position_dodge(.4)) + 
  scale_y_continuous(limits = c(-0.041,1.1))+
  ylab("A. monospora \n infection probability")+
  xlab("Food")+
  ggtitle("Experiment 2\nInfection Probability")+
  scale_x_discrete(labels=c("Ank","CYA"))+
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
    legend.position = c(0.7,0.7) # Adjusted argument
  ) +
  scale_fill_manual(
    name = "",
    labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
    values = c("#09753f", "red", "blue", "orange")
  ) +
  scale_shape_manual(
    name = "",
    labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
    values = c(21, 24, 21, 24)
  ) +
    geom_text(aes(x = 0.6, y = 0.95, label = "B)"), size = 6)
S4 # prints the plot

ggsave("metsch_infection_prob.png", plot = S4, width = 4, height = 4) # exports the plot

#### MATURE SPORES BURDEN ####
InfNoZero <- subset(subset_metsch_exposure, !is.na(infected) & mature > 0) #subests data to only include infected individuals
InfNoZero$parasite <- factor(InfNoZero$parasite) #parasite as factor
InfNoZero$food <- factor(InfNoZero$food) #food as factor

hist(InfNoZero$mature, #visualizes spore data with a histogram
     main = "Histogram of Metsch Spores",
     xlab = "Spores",
     ylab = "Frequency",)


#taking out the category wioth only one infected animal because you can't compare without error bars?
#making 3 separate cetegories by combining food and parasite and then removing the category with one infected made unique "treatment" with the food and parasite combinations

InfNoZero$treatment <- paste(InfNoZero$parasite,InfNoZero$food, sep="_") #creates new variable called treatment that is a combination of parasite and food
InfNoZero.analysis<-InfNoZero[!InfNoZero$treatment == "MicG/Metsch_CYA",] #removes MicG/Metsch_CYA because it only had one observation

mature_burden2<-glm(mature~treatment,data=InfNoZero.analysis) # creates generalized linear model w/ treatment
Anova(mature_burden2, type = 2) # prints stats

emmeans(mature_burden2,pairwise ~ treatment) #pairwise comparison

op <- par(mfrow = c(2, 2)) #diagnostic plots
plot(mature_burden2)             
par(op)


metsch_infection_matureSE<-summarySE(data = InfNoZero, measurevar = "mature", groupvars = c("parasite", "food"), na.rm=T) #summarizes data within treatments and calculates sd, se and ci


# Create the plot with all 4 groups for consistency but the analysis is only on treatment
S5 <- ggplot(metsch_infection_matureSE, aes(
  x = food, y = mature,
  fill = interaction(parasite, food),
  shape = interaction(parasite, food)
)) + 
  geom_errorbar(aes(ymin=mature-se, ymax=mature+se), colour="black", width=.3, position=position_dodge(0.4)) +
  geom_point(position=position_dodge(0.4), size=4)+
  scale_y_continuous(limits = c(0,200), expand = c(0,0))+
  ylab("A. monospora mature spore\nburden x10,000")+
  xlab("Food")+
  ggtitle("Experiemnt 2\nSpore Burden")+
  theme_classic()+
  theme(plot.margin = margin(20, 0, 0, 0),
        axis.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.title = element_text(size = 5, colour = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.position = c(0.7,0.8) 
  ) +
  scale_fill_manual(
    name = "",
    labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
    values = c("#09753f", "red", "blue", "orange")
  ) +
  scale_shape_manual(
    name = "",
    labels = c("AM Ank", "OP/AM Ank", "AM CYA", "OP/AM CYA"),
    values = c(21, 24, 21, 24)
  ) +
  geom_text(aes(x = 0.9, y = 185, label = "a"), size = 4) +
  geom_text(aes(x = 1.1, y = 90, label = "b"), size = 4) +
  geom_text(aes(x = 1.9, y = 65, label = "b"), size = 4) +
  geom_text(aes(x = 0.6, y = 180, label = "D)"), size = 6)

ggsave("spores.png", plot = S5, width = 4, height = 4) # exports the plot
S5

spore_plots_all <- grid.arrange(S2, S4, S3, S5, ncol = 2, nrow = 2)

ggsave("spore_plots_all.png", plot = spore_plots_all, width = 8, height = 8) # exports the plot

