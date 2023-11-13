#_____________________________________________31/05/2023___________________________________________________
# Here we use lipid concentrations to generate Tables 1 (lipid diversities and abundances), Figure S1 
# (Correlations between the abundances vs diversities of the various lipid subclasses in Day 1 and Day 
# 19 S06 males), Supplementary Tables S4, S5 and S6 and Figure S1 of Lipid MS1 
# lipids MS1 (Tables S4, S5 and S6 shows lipid species found in NLs,PLs, ELNLs, ELPLs and their abundances and 
# the number of replicates they were found in at both ages). We also calculate the PE/PC ratio at the end

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("gtools")
library("broom")
library("cowplot")
library("ggpubr")
library("gridExtra")
library("grid")
library("Rmisc")

rm(list=ls())

#_________________Outlier function to be used, when required______________________________ 
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#__________________Create a list of ambiguous lipid species from Shirleen's csv file______
ALS <- read_csv("data/Lipids to filter ploty.csv") %>% 
  dplyr::select(3) %>% 
  unique() %>% 
  c()
#______________ Load samples and remove outliers from batch02 ______________________________
# correct names of incorrectly labeled sample names
'%!in%' <- function(x,y)!('%in%'(x,y)) # a function to create opposite of %in% for filtering out

B01 <- read_csv("data/CD_results_Batch01_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))%>% 
  dplyr::filter(samples != "ct_n09")

B02 <- read_csv("data/CD_results_Batch02_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(samples != "S06_95") %>% 
  dplyr::mutate(samples = recode(samples,  syd_n87 = "syd_o87",syd_n104 = "syd_o104"))  

B03 <- read_csv("data/CD_results_Batch03_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>%
  dplyr::filter(samples != "syd_n134") %>% # remove outlier
  dplyr::mutate(samples = recode(samples, ct_n161 = "ct_o161", S06_171 = "s06_171",
                                 syd_160 = "syd_o160", syd_n138 = "syd_o138", syd_n169 = "syd_o169")) %>% 
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))

B04 <- read_csv("data/CD_results_Batch04_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(samples != "s06_217") %>% # remove outlier
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))

#_________________ calculate total area of all peaks in each sample and calculate their mean ___________________
F01 <- function(x){
  x %>% 
    dplyr::group_by(samples) %>% 
    dplyr::mutate(S.Area = sum(Area)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(samples, S.Area) %>% 
    unique() %>% 
    dplyr::mutate(M.S.Area = mean(S.Area)) %>% 
    dplyr::select(M.S.Area) %>% 
    slice(1)
}

MSA.B01 <- F01(x = B01)
MSA.B02 <- F01(x = B02)
MSA.B03 <- F01(x = B03)
MSA.B04 <- F01(x = B04)

NP.B01 <- B01 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B01))) %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

NP.B02 <- B02 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B02))) %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100)  

NP.B03 <- B03 %>% 
  dplyr::mutate(N.Area = Area)  %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

NP.B04 <- B04 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B04)))  %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

DF1 <- NP.B01 %>% 
  dplyr::bind_rows(NP.B02) %>% 
  dplyr::bind_rows(NP.B03) %>%
  dplyr::bind_rows(NP.B04) 
#______________________________________________________________________________________________________

#__________________________ Create main data frame _________________________________
Variables <- read_excel("data/Variables.xlsx") %>% 
  dplyr::select(ID, Weight, Time2) %>% 
  dplyr::rename(samples = ID) 

DF <- DF1 %>% 
  dplyr::left_join(., Variables, by = "samples") %>% 
  dplyr::filter(Weight < 7.5) %>% 
  dplyr::filter(Weight > 1) %>% 
  dplyr::select(Name, Weight, samples, N.Area, count, class, species, batch, DB_new, class_new) %>%
  dplyr::mutate(Bond = DB_new) %>% 
  dplyr::mutate(SubClass = class_new) %>% 
  dplyr::select(-DB_new, -class_new) %>% 
  tidyr::separate(samples, into = c("Line", "Cage"), remove = FALSE) %>% 
  dplyr::mutate(Line = replace(Line, Line == "S06", "s06")) %>% 
  dplyr::mutate(Class = class) %>% 
  dplyr::mutate(Batch = batch) %>%
  dplyr::select(-batch) %>% 
  dplyr::rename(Samples = samples) %>% 
  dplyr::mutate(Age = if_else(Batch %in% c("B01", "B02"), "19 Day", "1 Day")) %>% 
  dplyr::mutate(Time = str_extract(Cage, "[a-z]+")) %>% 
  dplyr::mutate(Time = replace_na(Time, "o")) %>% 
  dplyr::mutate(Time = if_else(Time %in% c("o"), "Old", "New")) %>% 
  dplyr::mutate(Time = replace(Time, Line == "s06", "Older")) %>% 
  dplyr::select(-count, -class) %>% 
  tidyr::separate(col = species, into = c("Carbon", NA), ":", remove = FALSE) %>% 
  dplyr::mutate(Carbon = as.double(Carbon))

#___________________

# Import standards data

df.old <- read_csv("Table_ConcAllStandards.csv") %>% # this is data from 1st run
  dplyr::select(Name, Class, area, Conc, Batch) %>%
  dplyr::filter(Class != "Cer" & Class != "FA" & Class != "SM" & Class != "PCp") %>% 
  dplyr::group_by(Name, Class, Conc, Batch) %>% 
  dplyr::summarise(Area = mean(area)) %>% # calculates total area for any class that had more than one standard 
  dplyr::ungroup() %>%
  unique()  %>% 
  dplyr::select(Conc, Name, Class, Area)

# Add manually calculate exp 2 CL areas to DF
DF.CL <- tibble(Conc = c(0.01, 0.10, 1.00, 10.00), 
                Name = rep("CL 72:4",4), Class = rep("CL", 4),
                Area = c(1023.129827,32570.02475,579994.8638,45993575.15))

DF.Stds <- rbind(df.old, DF.CL)

Regression <- DF.Stds %>% split(~Class) %>%
  map(~lm(log10(Area) ~ log10(Conc), data = .x)$coefficient) ### Obtain all the regression coefficient

Regression_coef <- tibble(Class2 = names(Regression) ,
                          bind_rows(lapply(Regression, as.data.frame.list))) %>% ### 
  dplyr::rename(Intercept =  "X.Intercept.",  Slope = "log10.Conc.") %>% 
  mutate(Class2 = if_else(Class2 =="PEp", "PE p", Class2))

#____________create main DF_____________________

DF.main <- DF %>%
  ungroup() %>% 
  mutate(Class2 = if_else(grepl(" e", SubClass), Class, SubClass)) %>% ## All ether lipids use ester lipid standards
  left_join(Regression_coef, by = "Class2") %>%
  mutate(Conc = 10^((log10(N.Area)-Intercept)/Slope)) %>%
  mutate(LogAreaCheck = Intercept + Slope*(log10(Conc)), AreaCheck = 10^LogAreaCheck) ### Just to double check

#_______1. Calculate percentage conc and number of lipid species detected for all lines___________

# Generate % sum peak area values & SE for Table 1 of manuscript

# complete fxn fills in zeros for samples that have missing lipids
Complete_DF <- DF.main %>% 
  complete(SubClass, nesting(Samples, Age, Line, Time), fill = list(Conc = 0)) %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Conc)

DF01 <- Complete_DF %>% 
  ungroup %>% 
  dplyr::group_by(Samples,SubClass) %>% 
  dplyr::mutate(Total.Conc.by.SAsC = sum(Conc)) %>% # by samples, age and subclass 
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>% 
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Total.Conc.by.SAsC, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.SAsC/Total.Conc.by.SA)*100) %>% 
  #run the codes till unique:total percentage should come to 100 for each individual sample to check: sum((DF01 %>% filter(Samples == "cbr_n04"))$Percentage)
  unique() %>% 
  dplyr::select(-Total.Conc.by.SAsC, -Total.Conc.by.SA, -Samples) %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Percentage") %>% 
  dplyr::select(Line, Time, Age, SubClass, Percentage, se)

# Check to see if % adds to 100 for an individual sample, run the above codes till unique() the total % should come to 100. 
# sum((DF01 %>% dplyr::filter(Samples == "s06_02"))$Percentage)

# Check to see if percentage adds to 100 for a Line, Age, Time
DF01 %>% 
  dplyr::filter(Line == "ct" & Time == "New" & Age == "19 Day") %>% 
  dplyr::select(Percentage) %>% 
  sum()

# format and save data for s06 only as csv
DF00 <- DF01 %>% 
  dplyr::mutate(Percentage = format(round(.$Percentage, 2), nsmall = 2)) %>% 
  dplyr::mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("PeSE", Percentage:se, sep = "±") %>% 
  pivot_wider(names_from = SubClass, values_from = PeSE) %>% 
  dplyr::filter(Line == "s06")

write.csv(DF00,"Lipid subclasses percent conc.csv")


# use DF to calculate % abundance & SE for each of the 4 lipid categories for s06 only for Table 1 of lipid MS1

# add another column to DF that has the following 4 lipid categories
phospholipids <- c("CL","PC", "PE", "PG", "PI", "PS","LPC")
neutrallipids <- c("DG", "TG")
PLetherlipids <- c("PC e", "PE e", "PE p", "PS e")
NLetherlipids <- c("DG e", "TG e")

DF_type <- Complete_DF %>%  
  ungroup %>% 
  dplyr::mutate(Type = ifelse(SubClass %in% phospholipids,"phospholipids", 
  ifelse(SubClass %in% PLetherlipids,"PLetherlipids", ifelse(SubClass %in% NLetherlipids,"NLetherlipids", "neutrallipids")))) %>% 
  dplyr::group_by(Samples, Type) %>% 
  dplyr::mutate(Total.Conc.by.Type = sum(Conc)) %>% # by samples, age and Type
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>% 
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-Conc) %>%
  dplyr::select(Samples, Age, Line, Time, Type, Total.Conc.by.Type, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.Type/Total.Conc.by.SA)*100) %>% # percentage should come to 100 for each individual sample to check: sum((DF_type %>% filter(Samples == "cbr_n04"))$Percentage)
  unique() %>% 
  dplyr::select(-Total.Conc.by.Type, -Total.Conc.by.SA) %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Percentage") %>% 
  dplyr::select(Line, Time, Age, Type, Percentage, se)

# Check to see if Percentage Type adds to 100 for a Line, Age, Time
DF_type %>% 
  dplyr::filter(Line == "s06" & Time == "Older" & Age == "1 Day") %>% 
  dplyr::select(Percentage) %>% 
  sum()

# format and save as csv
DF_type2 <- DF_type %>% 
  dplyr::mutate(Percentage = format(round(.$Percentage, 2), nsmall = 2)) %>% 
  dplyr::mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("PeSE", Percentage:se, sep = "±") %>%
  pivot_wider(names_from = Type, values_from = PeSE) %>% 
  dplyr::filter(Line == "s06")

write.csv(DF_type2,"Lipid catergories percent conc.csv")

#________2. Calculate percentage conc and number of lipid species for s06 to generate Tables S4, S5, S6_____

# First generate % conc values & SE for each of the lipid species present 
# NB: for the complete function used to fill in zeros we will have SubClass and Name inside the 1st nesting unlike above when calculating the % abundance of each SubClass, we only had SubClass in the 1st nesting

S06_lipids <- DF.main %>%
  ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Name, Conc) %>% 
  dplyr::filter(Line == "s06" & Samples != "s06_158")%>% # remove this outlier) %>% 
  complete(nesting(SubClass,Name), nesting(Samples, Age, Line, Time), fill = list(Conc = 0)) %>% 
  dplyr::group_by(Samples, SubClass, Name) %>% 
  dplyr::mutate(Total.Conc.by.SAsC = sum(Conc)) %>% # by samples, age and subclass, species 
  dplyr::ungroup() %>% 
  #dplyr::group_by(Samples,SubClass) %>% # sum area for a sample for each subclass
  dplyr::group_by(Samples) %>% # sum area for a sample regardless of subclass (we want this code grouping by samples so we find overall average abundance over 15/16 replicates not only in samples where a particular lipid species is present)
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Name, Total.Conc.by.SAsC, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.SAsC/Total.Conc.by.SA)*100) %>% # for each sample the lipid specie sum to 100% sum((S06_lipids %>% filter(Samples == "s06_02"))$Percentage) 
  unique() %>% 
  dplyr::group_by(Name, Age) %>% 
  dplyr::mutate(reps = sum(Percentage>0)) %>% 
  dplyr::select(-Total.Conc.by.SAsC, -Total.Conc.by.SA, -Samples) %>%
  Rmisc::summarySE(groupvars = c("Line", "Time", "Age", "SubClass", "Name", "reps"), measurevar = "Percentage") %>% # when we take mean we are summarising for whole s06 line and not for sample
  dplyr::select(Line, Time, Age, SubClass, Name, Percentage, se, N, reps) 
  
# check to see the sum of percentages within a subclass should be same as the one for that subclass percentage abundance obtained in Table 1 of lipids MS1
sum((S06_lipids %>% filter(Age == "19 Day" & SubClass == "PC"))$Percentage)

# format and save file to generate tables S4, S5, S6
S06_lipids01 <- S06_lipids %>% 
  dplyr::filter(Line == "s06") %>% 
  dplyr::mutate(Percentage = format(round(.$Percentage, 2), nsmall = 3)) %>%
  dplyr::select(Name, Age, SubClass,Percentage,reps) %>% 
  tidyr::unite("PeRep", Percentage:reps, sep = ",") %>%
  pivot_wider(names_from = Age, values_from = PeRep) %>% 
  tidyr::unite("Abundance", "1 Day":"19 Day", sep = ")(") %>%
  dplyr::mutate(Abundance = paste0(Abundance, ")")) %>% 
  dplyr::mutate(Abundance = paste0("(", Abundance)) #%>% 
  #dplyr::filter(SubClass %in% c("DG e", "TG e", "PC e", "PE e", "PE p", "PS e"))

write.csv(S06_lipids01,"Final_Scripts/S06_lipid species and abundance.csv") 

# check to see which lipids in S06 appears only at age 1 day and only at age 19 day I have included column N.Area
df_1day <- DF %>% ungroup()%>% dplyr::filter(Line == "s06" & Age == "1 Day" & Samples != "s06_158")%>% dplyr::select(Name, SubClass) 
df_19day <- DF %>% ungroup()%>% dplyr::filter(Line == "s06" & Age == "19 Day" & Samples != "s06_158")%>% dplyr::select(Name, SubClass)

# The R function setdiff is used to find elements which are in the first object but not in the second object

diff1 <- setdiff(df_1day, df_19day) #this prints out those lipid species that appear in 1 day but not in 19day
diff19 <- setdiff(df_19day, df_1day) #this prints out those lipid species that appear in 19day but not in 1day

# find lipids that are common at both ages
commonlipids <- intersect(df_1day, df_19day)

write.csv(diff1_SA,"diff1_SA.csv")  
write.csv(diff19_SA,"diff19_SA.csv")
write.csv(commonlipids_SA,"commonlipids_SA.csv") 

#_________

# we can also find the number of replicates in which each of the lipid species are found using codes below
rep <- DF %>% 
  dplyr::filter(Line == "s06") %>% 
  ungroup() %>% 
  group_by(Name, Age) %>% 
  dplyr::summarise(count = n())

# join the column count that contains info on nbr of samples that have a lipid species present

#diff1_rep <- inner_join(diff1_SA, rep) %>% dplyr::select (-se)# df for lipids found only at 1day
#diff19_rep <- inner_join(diff19_SA, rep) %>% dplyr::select (-se) # df for lipids found only at 1day
#commonlipids_rep <- inner_join(commonlipids_SA, rep) %>% dplyr::select (-se) # df for lipids found only at both ages

#______3. Generate Figure 2 for lipids MS1 showing correlations between the abundances vs diversities of lipid subclasses___________

# import data to be used for plotting
DF_diversity <- read_excel("S06_lipiddiversity.xlsx") %>%
  dplyr::select(-se)

# make correlation plot 

plot <- DF_diversity %>% 
  dplyr::mutate(Subclass = recode(Subclass, "PC e" = "PCe", "PE e" = "PEe", "PE p" = "PEp", "PS e" = "PSe", "DG e" = "DGe", "TG e" = "TGe")) %>% #remove space between th etherlipids
  dplyr::mutate(Subclass = factor(as.factor(Subclass), levels = c("CL","PC","PE","PG","PI","PS","LPC","DG","TG", "PCe","PEe","PEp","PSe", "DGe","TGe"))) %>%
  dplyr::mutate(Class = factor(as.factor(Class), levels = c("neutral lipids","phospholipids", "ether neutral lipids","ether phospholipids"))) %>%
  dplyr::mutate(Age = recode(Age, "1 Day" = "Day 1", "19 Day" = "Day 19")) %>% 
  ggplot(aes(log10(species), log10(Percentage.Conc), colour = Class))+
  geom_point()+ 
  geom_text_repel(aes(label = paste0(Subclass)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  #ggtitle("Correlation plot of lipid abundance versus lipid diversity")+
  labs(x= "Diversity (log10 number of compounds identified)", y = "Abundance (log10 % concentration)")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "bottom") + 
  #theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) +
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #facet_grid(~ Age)+
  facet_grid(Age ~ .)+
  theme(strip.background =element_rect(fill="Black"))+
  theme(strip.text = element_text(colour = 'white', size = 8))+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3) # cor() fxn computes the correlation coefficient
#NB:when we separate the plot by age there is less than 3 datapoints in one of the subclass so doesn't give r squared value

ggsave(plot = plot, width = 6, height = 4, units = "in", dpi = 300,filename = "FigureS1.jpg")

#_____________________________________________________________

# calculate PE/PC ratio for s06

PEPC <- DF.main %>%
  ungroup() %>% 
  dplyr::select(Samples, Line, Time, Age, SubClass, Conc) %>% 
  dplyr::filter(SubClass %in% c("PE", "PC")) %>% 
  dplyr::group_by(Samples, SubClass) %>%
  dplyr::mutate(S.Conc.C = sum(Conc)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-Conc) %>% 
  unique() %>% 
  tidyr::pivot_wider(values_from = S.Conc.C, names_from = SubClass) %>% 
  dplyr::mutate(Ratio1 = PE/PC) %>%
  dplyr::select(-"PE", -"PC") %>% 
  summarySE(groupvars = c("Line", "Time", "Age"), measurevar = "Ratio1") %>% 
  dplyr::select(Line, Time, Age, Ratio1, se) %>% 
  mutate(Ratio1 = format(round(.$Ratio1, 4), nsmall = 4)) %>% 
  mutate(se = format(round(.$se, 4), nsmall = 4)) %>% 
  tidyr::unite("PE/PC.SE", Ratio1:se, sep = "±")%>% 
  dplyr::filter(Line == "s06")

#_____________________________END___________________________

