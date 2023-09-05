#____________________________________________________________02/06/2023__________________________________________________________

# Here we generate bar plots to observe correlations in the distributions across side chain length, double bond and joint categories across 
# subclasses to create Panel A and separating acyl chains in ether lipids with "e" and "p" to generate Panel B of Figure 3 lipid MS1.
# The mean proportions (% abundances) of individual chains that were short vs medium vs long, and saturated, vs mono-unsaturated vs 
# polyunsaturated was generated from S05.R (all lipids) and S06.R (all lipids plus ether lipids with e and p separated)
# Removed PSe and DGe because the acyl chains could not be resolved, we just end up with a long single polyunsaturated acyl chain

library("tidyverse")
library("ggpubr")
library("readxl")
library("gtools")
library("cowplot")
library("gridExtra")
library("grid")
library("ggpattern") # since package ggpattern is not on cran we install it from github using remotes package #install.packages("remotes")#remotes::install_github("coolbutuseless/ggpattern")

rm(list=ls())

#_______________Prepare figures for lipids subclasses for Panel A_________

SC <- read_excel("SC_only.xlsx") %>% 
  dplyr::select(-Time, -LineTimeAge, -SE, -n) %>%
  dplyr::rename(Chainlength = AcylClass) %>% 
  dplyr::filter(SubClass %in% c("DG","TG","CL","PC","PE","PI","PG","PS","LPC"))

Bond <- read_excel("Bond_only.xlsx") %>% 
  dplyr::select(-Time, -LineTimeAge, -SE, -n) %>% 
  dplyr::rename(Bonds = AcylClass) %>% 
  dplyr::filter(SubClass %in% c("DG","TG","CL","PC","PE","PI","PG","PS","LPC"))

SC_Bond <- read_excel("SC_Bondcombined.xlsx") %>% 
  dplyr::select(-Time, -LineTimeAge, -SE, -n) %>% 
  dplyr::rename(SCBond = AcylClass)%>% 
  dplyr::filter(SubClass %in% c("DG","TG","CL","PC","PE","PI","PG","PS","LPC"))

#__________________

# Make barplots use stacks instead of individual columns as above)

# Plot for acyl chain lengths without legends
SCplot <- SC %>% 
  dplyr::mutate(SubClass = factor(as.factor(SubClass), levels = c("DG","TG","CL","PC","PE","PG","PI","PS","LPC"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>%
  dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass, MEAN, fill = Chainlength))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("green", "darkcyan", "red")) +
  labs(#title="Proportion of individual chains in each subclass that were short (S) vs medium (M) vs long (L)",
    x= "Lipid subclasses",                                   
    y="Proportion of acyl chain lengths")+                  
  #fill = "Acyl chain lengths")+
  theme_bw()+
  facet_wrap(~ Age) + 
  #theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_blank())+   
  theme(legend.position="none")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 10))+
  theme(legend.title = element_blank())+ 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))

ggsave(plot = SCplot, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/SCplot.jpg")


# Plot for acyl chain lengths with legends
SCplot1 <- SC %>% 
  dplyr::mutate(SubClass = factor(as.factor(SubClass), levels = c("DG","TG","CL","PC","PE","PG","PI","PS","LPC"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>%
  dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass, MEAN, fill = Chainlength))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("green", "darkcyan", "red")) +
  labs(title="Proportion of individual chains in each subclass that were short (S) vs medium (M) vs long (L)",
       x= "Lipid subclasses",                                   
       y="Proportion of acyl chain lengths",                  
       fill = "Acyl chain lengths")+
  theme_bw()+
  facet_wrap(~ Age) + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold"))+  
  #theme(legend.position="bottom")+ #remove legend title
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (10))) + theme(axis.title = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))

ggsave(plot = SCplot1, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/SCplot1.jpg")


# Plot for degrees of saturation without legends

Bondplot <- Bond %>% 
  dplyr::mutate(SubClass = factor(as.factor(SubClass), levels = c("DG","TG","CL","PC","PE","PG","PI","PS","LPC"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>%
  dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass, MEAN, fill = Bonds))+
  geom_col_pattern(
    #aes(pattern = Bonds, pattern_angle = Bonds, pattern_spacing = Bonds),
    aes(pattern = Bonds), 
    fill            = 'white',
    colour          = 'black', 
    #pattern_density = 0.45,   # fraction of the filled area which should be covered by the pattern, striping is increased to 50% of the fill area
    pattern_density = 0.05,
    pattern_spacing = 0.03,   # determines how far away individual elements are from each other.
    #pattern_fill    = 'black',
    pattern_colour  = 'black') +
  scale_pattern_manual(values=c('none', 'stripe', 'weave')) + # use no pattern for 0 DB, 'stripe’ pattern for 1 DB and ‘weave’ pattern for more than 1 DB
  theme_bw() +
  labs(#title="Proportion of individual chains in each subclass that were saturated (0), vs mono-unsaturated (1) vs polyunsaturated (X)",
    y="Proportion of double bonds")+           
  #fill = "Degrees of unsaturation")+
  theme_bw()+
  facet_wrap(~ Age) + 
  theme(legend.position="none")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + theme(axis.title = element_text(size = 10))+
  theme(legend.title = element_blank())+ 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))

ggsave(plot = Bondplot, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/Bondplot.jpg")

# Plot for degrees of saturation with legends

Bondplot1 <- Bond %>% 
  dplyr::mutate(SubClass = factor(as.factor(SubClass), levels = c("DG","TG","CL","PC","PE","PG","PI","PS","LPC","TG e","PC e","PE e","PE p"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>%
  #dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass, MEAN, fill = Bonds))+
  geom_col_pattern(
    #aes(pattern = Bonds, pattern_angle = Bonds, pattern_spacing = Bonds),
    aes(pattern = Bonds), 
    fill            = 'white',
    colour          = 'black', 
    #pattern_density = 0.45,   # fraction of the filled area which should be covered by the pattern, striping is increased to 50% of the fill area
    pattern_density = 0.05,
    pattern_spacing = 0.03,   # determines how far away individual elements are from each other.
    #pattern_fill    = 'black',
    pattern_colour  = 'black') +
  scale_pattern_manual(values=c('none', 'stripe', 'weave')) + # use no pattern for 0 DB, 'stripe’ pattern for 1 DB and ‘weave’ pattern for more than 1 DB
  theme_bw() +
  labs(title="Proportion of individual chains in each subclass that were saturated (0), vs mono-unsaturated (1) vs polyunsaturated (X)",
       x= "Lipid subclasses",                    
       y="Proportion of double bonds",          
       fill = "Degrees of unsaturation")+
  theme_bw()+
  facet_wrap(~ Age) + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"))+ 
  theme(legend.position="bottom")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + theme(axis.title = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))

ggsave(plot = Bondplot1, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "Bondplot1.jpg")

# Plot for acyl chain length and saturation together

SCBondplot <- SC_Bond %>% 
  dplyr::mutate(SubClass = factor(as.factor(SubClass), levels = c("DG","TG","CL","PC","PE","PG","PI","PS","LPC"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>%
  dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass, MEAN, fill = SCBond))+
  geom_col_pattern(aes(pattern = SCBond),
                   colour          = 'black', 
                   #pattern_density = 0.45,     # fraction of the filled area which should be covered by the pattern, striping is increased to 50% of the fill area
                   pattern_density = 0.05,
                   pattern_spacing = 0.03,     # determines how far away individual elements are from each other.
                   pattern_colour  = 'black') +
  scale_pattern_manual(values = c(L0 = "none", L1 = "stripe", LX = "weave", M0 = "none", M1 = "stripe", MX = "weave", S0 = "none", S1 = "stripe", SX = "weave")) +
  scale_fill_manual(values=c(L0 = "green", L1 ="green", LX ="green", M0 = "darkcyan", M1 ="darkcyan", MX ="darkcyan", S0 = "red", S1 ="red", SX ="red")) +
  labs(#title="Proportion of acyl chains that were short & saturated (S0), short & monounsaturated (S1), 
    #short & polyunsaturated (SX), medium & saturated (M0), medium & monounsaturated (M1), medium & polyunsaturated (MX),
    #long & saturated (L0), long & monounsaturated (L1), and long & polyunsaturated (LX)",
    x= "Lipid subclasses",                                            
    y="Proportion of acyl chain length and double bonds")+               
  #fill = "Acyl chain length and degrees of unsaturation")+
  theme_bw()+
  ##       fill = guide_legend(override.aes = list(pattern = "none")))+
  facet_wrap(~ Age) + 
  #theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_blank())+ 
  theme(legend.position="none")+ 
  #guides(fill = guide_legend(nrow = 1))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + theme(axis.title = element_text(size = 8))+
  theme(legend.title = element_blank())+ 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))

ggsave(plot = SCBondplot, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/SCBondplot.jpg")

# Plot for acyl chain length and saturation together with legends

SCBondplot1 <- SC_Bond %>% 
  dplyr::mutate(SubClass = factor(as.factor(SubClass), levels = c("DG","TG","CL","PC","PE","PG","PI","PS","LPC","TG e","PC e","PE e","PE p"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>%
  #dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass, MEAN, fill = SCBond))+
  geom_col_pattern(aes(pattern = SCBond),
                   colour          = 'black', 
                   #pattern_density = 0.45,     # fraction of the filled area which should be covered by the pattern, striping is increased to 50% of the fill area
                   pattern_density = 0.05,
                   pattern_spacing = 0.03,     # determines how far away individual elements are from each other.
                   pattern_colour  = 'black') +
  scale_pattern_manual(values = c(L0 = "none", L1 = "stripe", LX = "weave", M0 = "none", M1 = "stripe", MX = "weave", S0 = "none", S1 = "stripe", SX = "weave")) +
  scale_fill_manual(values=c(L0 = "green", L1 ="green", LX ="green", M0 = "darkcyan", M1 ="darkcyan", MX ="darkcyan", S0 = "red", S1 ="red", SX ="red")) +
  labs(title="Proportion of acyl chains that were short & saturated (S0), short & monounsaturated (S1), 
    short & polyunsaturated (SX), medium & saturated (M0), medium & monounsaturated (M1), medium & polyunsaturated (MX),
    long & saturated (L0), long & monounsaturated (L1), and long & polyunsaturated (LX)",
       x= "Lipid subclasses",                                            
       y="Proportion of acyl chain length and double bonds",               
       fill = "Acyl chain length and degrees of unsaturation")+
  theme_bw()+
  guides(pattern = guide_legend(override.aes = list(fill = "white")), #remove the pattern from the fill guide 
         fill = guide_legend(override.aes = list(pattern = "none")))+
  facet_wrap(~ Age) + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_blank())+ 
  theme(legend.position="bottom")+ 
  #guides(fill = guide_legend(nrow = 1))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (6))) + theme(axis.title = element_text(size = 6))+
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))

ggsave(plot = SCBondplot1, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/SCBondplot1.jpg")

#_____________

# arrange plots 

#p <- grid.arrange(SC.plot, Bond.plot, SCBond.plot, ncol=1,labels = "AUTO") #gridExtra package
#p1 <- ggarrange(SC.plot, Bond.plot, SCBond.plot, ncol=1, nrow =3,labels = "AUTO") #ggpubr package
#p2 <- plot_grid(SC.plot, Bond.plot, SCBond.plot, ncol=1, labels = "AUTO") #cowplot package 


#_______________Prepare figures for ether lipids separating ether and ester bonds for Panel B_________

rm(list=ls())

SCEL <- read_excel("SC_onlyEL.xlsx") %>%
  dplyr::mutate(Age2 = recode(Age2,"1 Day ester" = "Day 1 acyl", "19 Day ester" = "Day 19 acyl",
                              "1 Day ether" = "Day 1 alkyl/alkenyl", "19 Day ether" = "Day 19 alkyl/alkenyl")) %>%
  dplyr::select(-Time, -LineTimeAge, -SE, -n) %>%
  dplyr::rename(Chainlength = AcylClass) 

BondEL <- read_excel("Bond_onlyEL.xlsx") %>% 
  dplyr::mutate(Age2 = recode(Age2,"1 Day ester" = "Day 1 acyl", "19 Day ester" = "Day 19 acyl",
                              "1 Day ether" = "Day 1 alkyl/alkenyl", "19 Day ether" = "Day 19 alkyl/alkenyl")) %>%
  dplyr::select(-Time, -LineTimeAge, -SE, -n, -SubClass) %>% 
  dplyr::rename(Bonds = AcylClass)

SC_BondEL <- read_excel("SC_BondcombinedEL.xlsx") %>% 
  dplyr::mutate(Age2 = recode(Age2,"1 Day ester" = "Day 1 acyl", "19 Day ester" = "Day 19 acyl",
                              "1 Day ether" = "Day 1 alkyl/alkenyl", "19 Day ether" = "Day 19 alkyl/alkenyl")) %>%
  dplyr::select(-Time, -LineTimeAge, -SE, -n, -SubClass) %>% 
  dplyr::rename(SCBond = AcylClass)

#__________________

# Make barplots use stacks instead of individual columns as above)

# Plot for acyl chain lengths without legends
# SCEL %>% select(SubClass2) %>% unique()

SCplotEL <- SCEL %>% 
  dplyr::mutate(SubClass2 = factor(as.factor(SubClass2), levels = c("1 Day_TG e_ester", "1 Day_PC e_ester","1 Day_PE e_ester","1 Day_PE p_ester",
                                                                    "19 Day_TG e_ester","19 Day_PC e_ester","19 Day_PE e_ester", "19 Day_PE p_ester", 
                                                                    "1 Day_TG e_ether", "1 Day_PC e_ether","1 Day_PE e_ether","1 Day_PE p_ether",
                                                                    "19 Day_TG e_ether","19 Day_PC e_ether",  "19 Day_PE e_ether","19 Day_PE p_ether"))) %>%
  #dplyr::mutate(Age2 = recode(Age2, "1 Day ester" = "Day 1 ester", "19 Day ester" = "Day 19 ester",
  #                           "1 Day ether" = "Day 1 ether", "19 Day ether" = "Day 19 ether")) %>%
  dplyr::mutate(Age2 = factor(as.factor(Age2), levels = c("Day 1 alkyl/alkenyl","Day 19 alkyl/alkenyl","Day 1 acyl","Day 19 acyl"))) %>%
  #dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass2, MEAN, fill = Chainlength))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("green", "darkcyan", "red")) +
  labs(#title="Proportion of individual chains in each subclass that were short (S) vs medium (M) vs long (L)",
    x= "Lipid subclasses",                                   
    y="Proportion of acyl chain lengths")+                  
  theme_bw()+
  facet_grid(~ Age2, scales = "free", space = "free")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 10))+
  theme(legend.title = element_blank())+ 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))#+
  #theme(axis.text.x = element_blank())#+ theme(legend.position = "none")

ggsave(plot = SCplotEL, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/SCplotEL.jpg")

# Plot for degrees of saturation without legends

BondplotEL <- BondEL %>% 
  dplyr::mutate(SubClass2 = factor(as.factor(SubClass2), levels = c("1 Day_TG e_ester", "19 Day_TG e_ester", "1 Day_TG e_ether", "19 Day_TG e_ether",
                                                                    "1 Day_PC e_ester", "19 Day_PC e_ester", "1 Day_PC e_ether", "19 Day_PC e_ether", 
                                                                    "1 Day_PE e_ester", "19 Day_PE e_ester", "1 Day_PE e_ether", "19 Day_PE e_ether",
                                                                    "1 Day_PE p_ester", "19 Day_PE p_ester", "1 Day_PE p_ether", "19 Day_PE p_ether"))) %>%
  #dplyr::mutate(Age2 = recode(Age2, "1 Day ester" = "Day 1 ester", "19 Day ester" = "Day 19 ester",
  #                            "1 Day ether" = "Day 1 ether", "19 Day ether" = "Day 19 ether")) %>%
  dplyr::mutate(Age2 = factor(as.factor(Age2),  levels = c("Day 1 alkyl/alkenyl","Day 19 alkyl/alkenyl","Day 1 acyl","Day 19 acyl"))) %>%
  dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass2, MEAN, fill = Bonds))+
  geom_col_pattern(aes(pattern = Bonds), 
                   fill            = 'white',
                   colour          = 'black', 
                   #pattern_density = 0.45,   # fraction of the filled area which should be covered by the pattern, striping is increased to 50% of the fill area
                   pattern_density = 0.05,
                   pattern_spacing = 0.03,   # determines how far away individual elements are from each other.
                   #pattern_fill    = 'black',
                   pattern_colour  = 'black') +
  scale_pattern_manual(values=c('none', 'stripe', 'weave')) + # use no pattern for 0 DB, 'stripe’ pattern for 1 DB and ‘weave’ pattern for more than 1 DB
  theme_bw() +
  labs(#title="Proportion of individual chains in each subclass that were saturated (0), vs mono-unsaturated (1) vs polyunsaturated (X)",
    y="Proportion of double bonds")+           
  theme_bw()+
  facet_grid(~ Age2, scales = "free", space = "free")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + theme(axis.title = element_text(size = 10))+
  theme(legend.title = element_blank())+ 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))#+
  #theme(axis.text.x = element_blank())+ theme(legend.position = "none")

ggsave(plot = BondplotEL, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/BondplotEL.jpg")


# Plot for acyl chain length and saturation together

SCBondplotEL <- SC_BondEL %>% 
  dplyr::mutate(SubClass2 = factor(as.factor(SubClass2), levels = c("1 Day_TG e_ester", "19 Day_TG e_ester", "1 Day_TG e_ether", "19 Day_TG e_ether",
                                                                    "1 Day_PC e_ester", "19 Day_PC e_ester", "1 Day_PC e_ether", "19 Day_PC e_ether", 
                                                                    "1 Day_PE e_ester", "19 Day_PE e_ester", "1 Day_PE e_ether", "19 Day_PE e_ether",
                                                                    "1 Day_PE p_ester", "19 Day_PE p_ester", "1 Day_PE p_ether", "19 Day_PE p_ether"))) %>%
  # dplyr::mutate(Age2 = recode(Age2, "1 Day ester" = "Day 1 ester", "19 Day ester" = "Day 19 ester",
  #                             "1 Day ether" = "Day 1 ether", "19 Day ether" = "Day 19 ether")) %>%
  dplyr::mutate(Age2 = factor(as.factor(Age2), levels = c("Day 1 alkyl/alkenyl","Day 19 alkyl/alkenyl","Day 1 acyl","Day 19 acyl"))) %>%
  #dplyr::filter(Line == "s06") %>% 
  unique() %>%    
  ggplot(aes(SubClass2, MEAN, fill = SCBond))+
  geom_col_pattern(aes(pattern = SCBond),
                   colour          = 'black', 
                   #pattern_density = 0.45,     # fraction of the filled area which should be covered by the pattern, striping is increased to 50% of the fill area
                   pattern_density = 0.05,
                   pattern_spacing = 0.03,     # determines how far away individual elements are from each other.
                   pattern_colour  = 'black') +
  scale_pattern_manual(values = c(L0 = "none", L1 = "stripe", LX = "weave", M0 = "none", M1 = "stripe", MX = "weave", S0 = "none", S1 = "stripe", SX = "weave")) +
  scale_fill_manual(values=c(L0 = "green", L1 ="green", LX ="green", M0 = "darkcyan", M1 ="darkcyan", MX ="darkcyan", S0 = "red", S1 ="red", SX ="red")) +
  labs(#title="Proportion of acyl chains that were short & saturated (S0), short & monounsaturated (S1), 
    #short & polyunsaturated (SX), medium & saturated (M0), medium & monounsaturated (M1), medium & polyunsaturated (MX),
    #long & saturated (L0), long & monounsaturated (L1), and long & polyunsaturated (LX)",
    x= "Lipid subclasses",                                            
    y="Proportion of acyl chain length and double bonds")+               
  theme_bw()+
  facet_grid(~ Age2, scales = "free", space = "free")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + theme(axis.title = element_text(size = 8))+
  theme(legend.title = element_blank())+ 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.3))
  #theme(axis.text.x = element_blank())+ theme(legend.position = "none")

ggsave(plot = SCBondplotEL, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/SCBondplotEL.jpg")

#_____________________________________END_________________________
