#______________________________________22/05/2022____________________

# Here we make plots for all the standards used for the lipid quantification

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("cowplot")

rm(list=ls())

#_____________________________

# Make scatter plots for standards from 1st run

df.old <- read_csv("Table_ConcAllStandards.csv") %>% # this is data from 1st run
  dplyr::select(Name, Class, area, Conc, Batch) %>%
  dplyr::filter(Class != "Cer" & Class != "FA" & Class != "SM" & Class != "PCp") %>% 
  dplyr::group_by(Name, Class, Conc, Batch) %>% 
  dplyr::summarise(Area = mean(area)) %>% # calculates total area for any class that had more than one standard 
  dplyr::ungroup() %>%
  unique()  #this will give one datapoint for each standard used

#df.old %>% dplyr::select(Name) %>% unique()

plot <- df.old %>%
  dplyr::mutate(Class = factor(as.factor(Class), levels = c("DG","TG", "PC", "PE", "PG", "PI", "PS", "LPC", "PEp"))) %>%
  ggplot(aes(log10(Conc), log10(Area)))+
  geom_point(aes(colour = factor(Conc)))+
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  facet_grid(~Class)+
  labs(x= "Concentrations (log10 of concentration)", y = "Abundance (log10 peak area)")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size = 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 8)) +
  #theme(plot.title = element_text(size = 6)) +
  theme(strip.background =element_rect(fill="Black"))+
  theme(strip.text = element_text(colour = 'white', size = 8))+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 4) +
  stat_regline_equation(label.y = 8.5, geom = "label", size = 4)

ggsave(plot = plot, width = 12, height = 8, units = "in", dpi = 300,filename = "Final_ScriptsStandards/Figure1strun.jpg")

#_______________
df.new <- read_csv("Cardiolipin stds.csv") %>% # data from 2nd rum
  pivot_longer(cols = 2:6, names_to = "Name", values_to = "Area") %>% 
  tidyr::separate(col = Name, into = c("Class", NA), " ", remove = FALSE)%>% 
  dplyr::rename(Conc = std_conc)%>% 
  dplyr::mutate(Class = paste0(Class, "_2"))

df.old <-  df.old %>% #from above
  dplyr::mutate(Class = paste0(Class, "_1")) %>% 
  dplyr::select(Conc, Name, Class, Area)

DF <- rbind(df.new, df.old) %>% # merge the two sets of standard data
  dplyr::filter(Conc != "10" | Name != "TG 54:3" | Class != "TG_2") # remove the outlier in TG standard 2nd run

unique(DF$Class)

Combinedplot <- DF %>%
  dplyr::mutate(Class = factor(as.factor(Class), levels = c("DG_1","DG_2","TG_1","TG_2", "CL_2", "PC_1", "PC_2", "PE_1", "PE_2", "PG_1", "PI_1", "PS_1", "LPC_1", "PCp_1","PEp_1"))) %>%
  ggplot(aes(log10(Conc), log10(Area)))+
  geom_point(aes(colour = factor(Conc)))+
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  facet_grid(~Class)+
  labs(x= "Concentrations (log10 of concentration)", y = "Abundance (log10 peak area)")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size = 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 8)) +
  #theme(plot.title = element_text(size = 6)) +
  theme(strip.background =element_rect(fill="Black"))+
  theme(strip.text = element_text(colour = 'white', size = 8))+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 4) +
  stat_regline_equation(label.y = 8.5, geom = "label", size = 4)

ggsave(plot = Combinedplot, width = 17, height = 8, units = "in", dpi = 300,filename = "StandardsCombinedFigure.jpg")


# Add manually calculate exp 2 CL areas to DF
DF.CL <- tibble(Conc = c(0.01, 0.10, 1.00, 10.00), 
                Name = rep("CL 72:4",4), Class = rep("CL_1", 4),
                Area = c(1023.129827,32570.02475,579994.8638,45993575.15))
  

DF.CombinedExp <- rbind(DF, DF.CL)

Combinedplot2 <- DF.CombinedExp %>%
  dplyr::mutate(Class = factor(as.factor(Class), levels = c("DG_1","DG_2","TG_1","TG_2", "CL_1", "CL_2", "PC_1", "PC_2", "PE_1", "PE_2", "PG_1", "PI_1", "PS_1", "LPC_1", "PCp_1","PEp_1"))) %>%
  ggplot(aes(log10(Conc), log10(Area)))+
  geom_point(aes(colour = factor(Conc)))+
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  facet_grid(~Class)+
  labs(x= "Concentrations (log10 of concentration)", y = "Abundance (log10 peak area)")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size = 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 8)) +
  #theme(plot.title = element_text(size = 6)) +
  theme(strip.background =element_rect(fill="Black"))+
  theme(strip.text = element_text(colour = 'white', size = 8))+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 4) +
  stat_regline_equation(label.y = 8.5, geom = "label", size = 4)


#___________________________________________________________________________________________________
#___________________________________________________________________________________________________

#______Make the plots for "DG","TG", "CL","PC", "PE" with 1st (4 batches) and 2nd standard runs___________
df.new01 <- read_csv("Cardiolipin stds.csv") %>% # data from 2nd rum
  pivot_longer(cols = 2:6, names_to = "Name", values_to = "Area") %>% 
  tidyr::separate(col = Name, into = c("Class", NA), " ", remove = FALSE)%>% 
  dplyr::rename(Conc = std_conc)%>% 
  dplyr::mutate(Batch = 5) %>% 
  dplyr::select(Name, Class, Area, Conc, Batch) %>% 
  dplyr::filter(Conc != "10" | Name != "TG 54:3" | Class != "TG") # remove the outlier in TG standard 2nd run 

df.old01 <- read_csv("Table_ConcAllStandards.csv") %>% # this is data from 1st run
  dplyr::select(Name, Class, area, Conc, Batch) %>% 
  dplyr::group_by(Class,Conc, Batch) %>% 
  dplyr::mutate(Area = mean(area)) %>% # calculates overall means from the four batches and technical replicates
  dplyr::ungroup() %>% 
  dplyr::filter(Class != "Cer" & Class != "FA" & Class != "SM" & Class != "PCp") %>% 
  dplyr::select(Name, Class, Area, Conc, Batch) 

DF.Combined <- rbind(df.new01,df.old01) %>% 
  dplyr::filter(Class %in% c("DG" ,"TG" ,"PE","PC")) %>% 
  unique()

Batchplot <- DF.Combined %>%
  dplyr::mutate(Class = factor(as.factor(Class), levels = c("DG","TG", "PC", "PE"))) %>%
  ggplot(aes(log10(Conc), log10(Area), shape = Batch))+
  geom_point(aes(colour = factor(Conc)))+
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  facet_grid(~Class)+
  labs(x= "Concentrations (log10 of concentration)", y = "Abundance (log10 peak area)")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "bottom") + 
  theme(legend.title=element_text(size = 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 8)) +
  #theme(plot.title = element_text(size = 6)) +
  theme(strip.background =element_rect(fill="Black"))+
  theme(strip.text = element_text(colour = 'white', size = 8))#+
#stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 4) +
#stat_regline_equation(label.y = 8.5, geom = "label", size = 4)

ggsave(plot = Batchplot, width = 10, height = 8, units = "in", dpi = 300,filename = "Final_Scripts/Figurewithstandardruninbatches.jpg")

#_________________END_______________

